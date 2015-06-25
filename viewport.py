#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from os.path import dirname
from mol_f import make_iso_surf
import numpy as np

from PyQt4.QtCore import Qt
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
from OpenGL.GL import *
from OpenGL.GL.shaders import *
from OpenGL.arrays.vbo import *

class ViewPort(QGLWidget):

        ##################################################
        # CALLED UPON INITIALISATION
        ##################################################
        def __init__(self,parent):
                #init with antialiasing, transparency and OGLv3.3 core profile
                form=QGLFormat(QGL.SampleBuffers|QGL.AlphaChannel)
                form.setProfile(QGLFormat.CoreProfile)
                super(ViewPort,self).__init__(form)
                self.parent = parent
                self.config=self.parent.controller._config
                self.perspective=self.config['Proj-Perspective']
                self.showBonds = True
                self.showCell = True
                self.showPlane = False
                self.showSurf = False
                self.mouseMode = "Camera"
                self.AA = True
                self.xsh = 0
                self.ysh = 0
                self.rMatrix = QMatrix4x4()
                self.distance = 25
                self.initActions()

        def initActions(self):
                changeProj = QAction('Change &Projection',self)
                changeProj.setShortcut('p')
                changeProj.triggered.connect(self.changeState)
                toggleCell = QAction('Toggle &Cell',self)
                toggleCell.setShortcut('c')
                toggleCell.triggered.connect(self.changeState)
                antiAlias = QAction('Toggle &Antialiasing',self)
                antiAlias.setShortcut('a')
                antiAlias.triggered.connect(self.changeState)
                toggleBonds = QAction('Toggle &Bonds',self)
                toggleBonds.setShortcut('b')
                toggleBonds.triggered.connect(self.changeState)
                self.undoAction = QAction('Undo',self)
                self.undoAction.setShortcut('Ctrl+Z')
                self.undoAction.triggered.connect(self.undo)
                vMenu = self.parent.parent.menuBar().addMenu('&View')
                vMenu.addAction(changeProj)
                vMenu.addAction(toggleCell)
                vMenu.addAction(toggleBonds)
                vMenu.addAction(antiAlias)
                vMenu.addSeparator()
                vMenu.addAction(self.undoAction)

        def initializeGL(self):
                #Bind VAO (necessary for modern OpenGL)
                glBindVertexArray(glGenVertexArrays(1))
                #render only visible vertices
                glEnable(GL_DEPTH_TEST)
                #backface culling: render only front of vertex
                glEnable(GL_CULL_FACE)
                #enable transparency for selection
                glEnable(GL_BLEND)
                glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)

                #set line width for cell
                glLineWidth(2)
                glPointSize(2)

                #set background color
                self.qglClearColor(QColor(255,255,255,0))

                #find supported version and modify shaders accordingly
                glVersion=float(glGetString(GL_VERSION)[0:3])
                self.instanced = glVersion>=3.3
                if self.instanced:
                    self.glslv='#version 330\n'
                else:
                    self.glslv='#version 130\n'

                #add shaders:
                self.sphereShader=self.makeShader('vertexSpheres.vsh','fragmentSpheres.fsh')
                self.bondShader=self.makeShader('vertexBonds.vsh','fragmentBonds.fsh')
                self.lineShader=self.makeShader('vertexLines.vsh','fragmentLines.fsh')
                self.selectShader=self.makeShader('vertexSelect.vsh','fragmentSelect.fsh')
                self.planeShader=self.makeShader('vertexPlane.vsh','fragmentPlane.fsh')
                self.surfShader=self.makeShader('vertexSurf.vsh','fragmentSpheres.fsh')

                # load sphere
                sf=open(dirname(__file__)+'/sphere_model','r')
                self.sphereVBO = VBO(np.array(sf.readline().split(),'f'))
                sf.close()
                # load torus
                tf=open(dirname(__file__)+'/bond_model','r')
                self.torusVBO=VBO(np.array(tf.readline().split(),'f'))
                tf.close()

        def makeShader(self,vf,ff):
            s = QGLShaderProgram()
            with open(dirname(__file__)+'/'+vf) as f:
                v=self.glslv+f.read()
                s.addShaderFromSourceCode(QGLShader.Vertex,v)
            with open(dirname(__file__)+'/'+ff) as f:
                f=self.glslv+f.read()
                s.addShaderFromSourceCode(QGLShader.Fragment,f)
            return s

        ##########################################
        # Modify state of Visualization
        ##########################################

        def changeState(self):
                reason = self.sender().text()
                if reason=='Change &Projection':
                    self.perspective = not self.perspective
                    self.updateGL()
                elif reason=='Toggle &Cell':
                    self.showCell = not self.showCell
                    self.updateGL()
                elif reason=='Toggle &Bonds':
                    self.showBonds = not self.showBonds
                    self.updateGL()
                elif reason=='Toggle &Antialiasing':
                    if self.AA:
                            glDisable(GL_MULTISAMPLE)
                    elif not self.AA:
                            glEnable(GL_MULTISAMPLE)
                    self.AA = not self.AA
                    self.updateGL()

        def alignView(self):
            if self.sender().text()=='x':
                self.rMatrix = QMatrix4x4([0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1])
            elif self.sender().text()=='y':
                self.rMatrix = QMatrix4x4([1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1])
            elif self.sender().text()=='z':
                self.rMatrix.setToIdentity()
            self.updateGL()

        def setMouseMode(self,but):
            t=but.text()
            self.mouseMode = t
            if t == "Camera":
                self.setCursor(Qt.ArrowCursor)
            elif t == "Select":
                self.setCursor(Qt.CrossCursor)
            elif t == "Modify":
                self.setCursor(Qt.OpenHandCursor)

        #########################
        # Update stuff that's going to be drawn
        #########################

        def setMol(self,mol,mult):
                #prepare atoms and bonds for drawing
                if not mol:return
                #deactivate plane and volume data when
                #molecule changes
                if hasattr(self,'mol') and (self.mol is not mol):
                    self.showPlane = False
                    self.showSurf = False
                    if hasattr(self,'planeVBO'):
                        del self.planeVBO
                    if hasattr(self,'surfVBO'):
                        del self.surfVBO

                #check if molecule has undo-stack:
                undo = mol.get_undo()
                if undo:
                    self.undoAction.setText('Undo '+undo)
                    self.undoAction.setEnabled(True)
                else:
                    self.undoAction.setText('Undo')
                    self.undoAction.setDisabled(True)

                #clear old selection if present:
                if hasattr(self,'selVBO'):
                    del self.selVBO

                #save for interaction
                self.mol=mol
                self.mult=mult
                #local variables for convenience
                atoms = mol.get_all_atoms()
                pse = mol.pse
                vec = mol.get_vec()*mol.get_celldm()
                center = mol.get_center(self.config['Rotate-COM'])
                bonds = mol.get_bonds(self.config['Bond-Cut-Fac'])
                sel = mol.get_selection()

                #get bonds and calculate offsets
                if mult == [1,1,1]:
                        #only one (==no) offset
                        off = [-center]
                        edge=[7]
                else:
                        tmult = [1,1,1]
                        #save the multiplicators:
                        for i,j in enumerate(mult):
                                if j%2 == 0:
                                        tmult[i]=[x+0.5-j/2 for x in range(j)]
                                else:
                                        tmult[i]=[x-np.floor(j/2) for x in range(j)]
                        #generate offsets:
                        off=[(i*vec[0]+j*vec[1]+k*vec[2])-center for i in tmult[0]
                                for j in tmult[1] for k in tmult[2]]
                        #save binary representation of disabled pbc-bonds (b'zyx')
                        edge=[ (i==mult[0]) + ((j==mult[1])<<1) + ((k==mult[2])<<2)
                                for i in range(1,mult[0]+1)
                                for j in range(1,mult[1]+1)
                                for k in range(1,mult[2]+1)]

                #prepare bond VBOs
                self.bondPosVBO=[]
                #binary representation of enabled pbc directions (b'zyx')
                mult=np.sign(mult[0]-1)+np.sign(mult[1]-1)*2+np.sign(mult[2]-1)*4
                r=self.config['Bond-Radius']
                for j in [0,1,2,3,4,5,6,7]:
                        #check if pbc part is necessary
                        if mult&j!=j: continue
                        for i in bonds[j]:
                                #get positions of atoms
                                a = atoms[i[0]][1]+i[2][0]
                                b = atoms[i[1]][1]+i[2][1]
                                #save colors
                                c1 = pse[atoms[i[0]][0]][5:]
                                c2 = pse[atoms[i[1]][0]][5:]
                                #length-scaling-factor:
                                l = i[3]
                                #position of bond
                                pos= (a+b)/2
                                #rotate bond from x-axis d to bond-axis c
                                c = (a-b)/np.linalg.norm(a-b)
                                d = np.array([1,0,0],'f')
                                #check if parallel to x-axis
                                if np.all(np.equal(abs(c),d)):
                                        ax=np.array([0,1,0],'f')
                                        c=c[0]
                                        s=0
                                        ic=1-c
                                else:
                                        theta=np.arccos(np.dot(c,d))
                                        ax = -np.cross(c,d)
                                        ax=ax/np.linalg.norm(ax)
                                        #construct rotation matrix
                                        c=np.float(np.cos(theta))
                                        s=np.float(-np.sin(theta))
                                        ic=np.float(1.-c)
                                for idx,k in enumerate(off):
                                        if j>0 and edge[idx]&j!=0:
                                                continue
                                        self.bondPosVBO.append([l*(ic*ax[0]*ax[0]+c),l*(ic*ax[0]*ax[1]-s*ax[2]),l*(ic*ax[0]*ax[2]+s*ax[1]),0.,
                                                    r*(ic*ax[0]*ax[1]+s*ax[2]),r*(ic*ax[1]*ax[1]+c),r*(ic*ax[1]*ax[2]-s*ax[0]),0.,
                                                    r*(ic*ax[0]*ax[2]-s*ax[1]),r*(ic*ax[1]*ax[2]+s*ax[0]),r*(ic*ax[2]*ax[2]+c),0.,
                                                    pos[0]+k[0],pos[1]+k[1],pos[2]+k[2],1.,
                                                    c1[0],c1[1],c1[2],c1[3],
                                                    c2[0],c2[1],c2[2],c2[3]])
                if self.instanced:
                    self.bondPosVBO=VBO(np.array(self.bondPosVBO,'f'))

                #save atoms in VBOs
                self.atomsVBO=[(at[1]+j).tolist()+[pse[at[0]][3+self.config['Atom-Radius-VdW']],pse[at[0]][5],pse[at[0]][6],pse[at[0]][7],pse[at[0]][8]] for at in atoms for j in off]
                if self.instanced:
                    self.atomsVBO=VBO(np.array(self.atomsVBO,'f'))
                #check for selected atoms inside mult-range
                if sel:
                    self.selVBO=[]
                    for i in sel:
                        if all(i[1]<np.array(self.mult)):
                            at=atoms[i[0]]
                            pos = (at[1]+np.dot(i[1],vec)+off[0]).tolist()
                            self.selVBO.append(pos+[pse[at[0]][3+self.config['Atom-Radius-VdW']]*1.5,0.4,0.4,0.5,0.5])
                    if self.selVBO and self.instanced:
                        self.selVBO=VBO(np.array(self.selVBO,'f'))
                #make cell:
                null=np.zeros(3)
                celltmp=[null,vec[0],null,vec[1],null,vec[2],
                        vec[0],vec[0]+vec[1],vec[0],vec[0]+vec[2],
                        vec[1],vec[1]+vec[0],vec[1],vec[1]+vec[2],
                        vec[2],vec[2]+vec[0],vec[2],vec[2]+vec[1],
                        vec[0]+vec[1],vec[0]+vec[1]+vec[2],
                        vec[0]+vec[2],vec[0]+vec[1]+vec[2],
                        vec[1]+vec[2],vec[0]+vec[1]+vec[2]]
                self.cellVBO=VBO(np.array([i+j for j in off for i in celltmp],'f'))

                # save offset for plane and volume
                if self.instanced:
                    self.offVBO=VBO(np.array(off,'f'))
                else:
                    self.offVBO=off

                self.updateGL()

        ################################################
        # CREATE AND MANAGE PLANES AND SURFACES
        ################################################
        def toggleSurf(self):
            self.showSurf = not self.showSurf
            self.updateGL()

        def setSurf(self,sval,both=False):
            if sval==0:
                return
            vertices,nv = make_iso_surf(self.mol.get_vol(),sval,self.mol.get_vol_gradient())
            if both:
                vertices2,nv2=make_iso_surf(self.mol.get_vol(),-sval,self.mol.get_vol_gradient())
                self.surfVBO = VBO(np.concatenate((vertices[:,:,:nv].flatten('F'),vertices2[:,:,:nv2].flatten('F')),axis=0))
            else:
                self.surfVBO = VBO(vertices[:,:,:nv].flatten('F'))
            if self.showSurf:
                self.updateGL()

        def togglePlane(self):
                self.showPlane = not self.showPlane
                self.updateGL()

        def setPlane(self,ptype,pval):
                #volume data:
                #'x/y/z',int
                if ptype in 'xyz':
                    v = self.mol.get_vol()
                    vmin = v.min()
                    vdiff = v.max()-vmin
                    if ptype=='x':
                        pdat=v[pval,:,:]
                        pval = pval/v.shape[0]
                        p=[[pval,0,0],[pval,1,0],[pval,0,1],[pval,1,1]]
                    elif ptype=='y':
                        pdat=v[:,pval,:]
                        pval = pval/v.shape[1]
                        p=[[0,pval,0],[1,pval,0],[0,pval,1],[1,pval,1]]
                    elif ptype=='z':
                        pdat=v[:,:,pval]
                        pval = pval/v.shape[2]
                        p=[[0,0,pval],[1,0,pval],[0,1,pval],[1,1,pval]]
                    self.planeTex=np.array(map(lambda x:(x-vmin)/vdiff,pdat),'f')

                #crystal data:
                #'c',[tuple]
                elif ptype == 'c':
                    self.planeTex=np.array([[1.]],'f')
                    #catch undefined case
                    if pval.count(0) == 3:
                        if hasattr(self,'planeVBO'):
                            del self.planeVBO
                        self.updateGL()
                        return
                    elif pval[0] == 0:
                        if pval[1] == 0:
                            p=[[0,0,pval[2]],[1,0,pval[2]],[0,1,pval[2]],[1,1,pval[2]]]
                        elif pval[2] == 0:
                            p=[[0,pval[1],0],[1,pval[1],0],[0,pval[1],1],[1,pval[1],1]]
                        else:
                            p=[[0,0,pval[2]],[0,pval[1],0],[1,0,pval[2]],[1,pval[1],0]]
                    else:
                        if pval[1] == 0:
                            if pval[2]==0:
                                p=[[pval[0],0,0],[pval[0],1,0],[pval[0],0,1],[pval[0],1,1]]
                            else:
                                p=[[pval[0],0,0],[pval[0],1,0],[0,0,pval[2]],[0,1,pval[2]]]
                        elif pval[2] == 0:
                            p=[[pval[0],0,0],[pval[0],0,1],[0,pval[1],0],[0,pval[1],1]]
                        else:
                            p=[[pval[0],0,0],[0,pval[1],0],[0,0,pval[2]]]
                p=np.array(p,'f')
                #take care of negative hkl-values
                if ptype == 'c':
                    for i in range(3):
                        if pval[i]<0:
                            p[:,i]+=1

                #generate planeVBO
                vec=self.mol.get_vec()*self.mol.get_celldm()
                UV = [[0,0],[0,1],[1,0],[1,1]]
                self.planeVBO=VBO(np.array([np.dot(p[i],vec).tolist()+UV[i] for i in range(len(p))],'f'))

                if self.showPlane:
                    self.updateGL()
                return

        ################################################
        # CALLED UPON WINDOW RESIZE
        ################################################
        def resizeGL(self,width,height):
                #prevent divide by zero
                if height == 0: height = 1

                aspect = float(width)/float(height)
                #set projection matrix
                self.pMatrix = QMatrix4x4()
                self.pMatrix.setToIdentity()
                self.pMatrix.perspective(60.0,aspect,0.001,1000)
                #set orthogonal matrix:
                self.oMatrix = QMatrix4x4()
                self.oMatrix.setToIdentity()
                self.oMatrix.ortho(-10*aspect,10*aspect,-10,10,0.001,1000)

                #set viewport
                glViewport(0,0,width,height)

        ###############################################
        # CALLED UPON WINDOW UPDATE EVENT
        ##############################################

        def paintGL(self,select=False):
                #clear depth and color buffer:
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

                if not hasattr(self,'atomsVBO'): return

                #construct viewMatrix
                self.vMatrix = QMatrix4x4()
                self.vMatrix.lookAt(QVector3D(0,0,self.distance),QVector3D(0,0,0),QVector3D(0,1,0))
                self.vMatrix.translate(self.xsh,self.ysh,0)

                #TODO: orthogonal zooming needs fix
                #check for projection:
                if self.perspective:
                        self.proj = self.pMatrix
                else:
                        self.proj = self.oMatrix
                        #scale based on distance for zoom effect
                        self.vMatrix.scale(10./self.distance)
                #rendering:
                if select:
                        self.drawAtomsSelect()
                else:
                        self.drawAtoms()
                        if self.showBonds:
                                self.drawBonds()
                        if self.showCell:
                                self.drawCell()
                        if hasattr(self,'selVBO'):
                                self.drawSelection()
                        if self.showPlane and hasattr(self,'planeVBO'):
                                self.drawPlane()
                        if self.showSurf and hasattr(self,'surfVBO'):
                            self.drawSurf()

        def drawAtoms(self):
            self.sphereShader.bind()

            self.sphereVBO.bind()
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
            self.sphereVBO.unbind()

            self.sphereShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
            self.sphereShader.setUniformValue('rMatrix',self.rMatrix)
            self.sphereShader.setUniformValue('atom_fac',self.config['Atom-Factor'])

            if self.instanced:
                self.atomsVBO.bind()
                glEnableVertexAttribArray(1)
                glEnableVertexAttribArray(2)
                glEnableVertexAttribArray(3)
                glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
                glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.atomsVBO+12)
                glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,32,self.atomsVBO+16)
                glVertexAttribDivisor(1,1)
                glVertexAttribDivisor(2,1)
                glVertexAttribDivisor(3,1)
                self.atomsVBO.unbind()
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO)/3,len(self.atomsVBO))
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glDisableVertexAttribArray(3)
                glVertexAttribDivisor(1,0)
                glVertexAttribDivisor(2,0)
                glVertexAttribDivisor(3,0)
            else:
                for i in self.atomsVBO:
                    self.sphereShader.setUniformValue('position_modelspace',*i[0:3])
                    self.sphereShader.setUniformValue('scale_modelspace',i[3])
                    self.sphereShader.setUniformValue('color_input',*i[4:])
                    glDrawArrays(GL_TRIANGLES,0,len(self.sphereVBO)/3)

            #reset
            glDisableVertexAttribArray(0)
            self.sphereShader.release()

        def drawSurf(self):
            self.surfShader.bind()

            #send vertices
            self.surfVBO.bind()
            glEnableVertexAttribArray(0)
            glEnableVertexAttribArray(1)
            glEnableVertexAttribArray(2)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,36,None)
            glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,36,self.surfVBO+12)
            glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,36,self.surfVBO+24)
            self.surfVBO.unbind()

            self.surfShader.setUniformValue('volOff',*self.mol.get_vol_offset().tolist())
            self.surfShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
            self.surfShader.setUniformValue('cellVec',QMatrix3x3((self.mol.get_vec()*self.mol.get_celldm()).flatten()))
            self.surfShader.setUniformValue('rMatrix',self.rMatrix)

            glDisable(GL_CULL_FACE)

            if self.instanced:
                self.offVBO.bind()
                glEnableVertexAttribArray(3)
                glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,0,None)
                glVertexAttribDivisor(3,1)
                self.offVBO.unbind()
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.surfVBO)/9,len(self.offVBO))
                glVertexAttribDivisor(3,0)
                glDisableVertexAttribArray(3)
            else:
                for i in self.offVBO:
                    self.surfShader.setUniformValue('offset',*i)
                    glDrawArrays(GL_TRIANGLES,0,len(self.surfVBO)/9)

            glEnable(GL_CULL_FACE)
            glDisableVertexAttribArray(0)
            glDisableVertexAttribArray(1)
            glDisableVertexAttribArray(2)
            self.surfShader.release()

        def drawPlane(self):
            self.planeShader.bind()

            self.planeVBO.bind()
            glEnableVertexAttribArray(0)
            glEnableVertexAttribArray(1)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,20,None)
            glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,20,self.planeVBO+12)
            self.planeVBO.unbind()

            self.planeShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)

            #send texture
            ID = glGenTextures(1)
            glActiveTexture(GL_TEXTURE0)
            glBindTexture(GL_TEXTURE_2D,ID)
            glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
            glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
            glTexImage2D(GL_TEXTURE_2D,0,GL_RED,self.planeTex.shape[1],self.planeTex.shape[0],0,GL_RED,GL_FLOAT,self.planeTex)
            self.planeShader.setUniformValue('texSampler',0)

            glDisable(GL_CULL_FACE)

            if self.instanced:
                self.offVBO.bind()
                glEnableVertexAttribArray(2)
                glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,0,None)
                glVertexAttribDivisor(2,1)
                self.offVBO.unbind()
                glDrawArraysInstanced(GL_TRIANGLE_STRIP,0,len(self.planeVBO),len(self.offVBO))
                glVertexAttribDivisor(2,0)
                glDisableVertexAttribArray(2)
            else:
                for i in self.offVBO:
                    self.planeShader.setUniformValue('offset',*i)
                    glDrawArrays(GL_TRIANGLE_STRIP,0,len(self.planeVBO))

            glEnable(GL_CULL_FACE)
            glDeleteTextures(1)
            glDisableVertexAttribArray(0)
            glDisableVertexAttribArray(1)
            self.planeShader.release()

        def drawAtomsSelect(self):
            self.selectShader.bind()

            self.sphereVBO.bind()
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
            self.sphereVBO.unbind()

            self.selectShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
            self.selectShader.setUniformValue('atom_fac',self.config['Atom-Factor'])

            if self.instanced:
                self.atomsVBO.bind()
                glEnableVertexAttribArray(1)
                glEnableVertexAttribArray(2)
                glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
                glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.atomsVBO+12)
                glVertexAttribDivisor(1,1)
                glVertexAttribDivisor(2,1)
                self.atomsVBO.unbind()
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO)/3,len(self.atomsVBO))
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glVertexAttribDivisor(1,0)
                glVertexAttribDivisor(2,0)
            else:
                for j,i in enumerate(self.atomsVBO):
                    self.selectShader.setUniformValue('position_modelspace',*i[0:3])
                    self.selectShader.setUniformValue('scale_modelspace',i[3])
                    self.selectShader.setUniformValue('in_color',(j&0xFF)/255.,((j&0xFF00)>>8)/255.,((j&0xFF0000)>>16)/255.,1)
                    glDrawArrays(GL_TRIANGLES,0,len(self.sphereVBO)/3)

            glDisableVertexAttribArray(0)
            self.selectShader.release()

        def drawSelection(self):
            self.sphereShader.bind()

            self.sphereVBO.bind()
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
            self.sphereVBO.unbind()

            self.sphereShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
            self.sphereShader.setUniformValue('rMatrix',self.rMatrix)
            self.sphereShader.setUniformValue('atom_fac',self.config['Atom-Factor'])

            if self.instanced:
                self.selVBO.bind()
                glEnableVertexAttribArray(1)
                glEnableVertexAttribArray(2)
                glEnableVertexAttribArray(3)
                glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
                glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.selVBO+12)
                glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,32,self.selVBO+16)
                glVertexAttribDivisor(1,1)
                glVertexAttribDivisor(2,1)
                glVertexAttribDivisor(3,1)
                self.selVBO.unbind()
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO)/3,len(self.selVBO))
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glDisableVertexAttribArray(3)
                glVertexAttribDivisor(1,0)
                glVertexAttribDivisor(2,0)
                glVertexAttribDivisor(3,0)
            else:
                for i in self.selVBO:
                    self.sphereShader.setUniformValue('position_modelspace',*i[0:3])
                    self.sphereShader.setUniformValue('scale_modelspace',i[3])
                    self.sphereShader.setUniformValue('color_input',*i[4:])
                    glDrawArrays(GL_TRIANGLES,0,len(self.sphereVBO)/3)

            glDisableVertexAttribArray(0)
            self.sphereShader.release()

        def drawBonds(self):
            self.bondShader.bind()

            self.torusVBO.bind()
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
            self.torusVBO.unbind()

            self.bondShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
            self.bondShader.setUniformValue('rMatrix',self.rMatrix)

            if self.instanced:
                self.bondPosVBO.bind()
                for i in range(1,7):
                    glEnableVertexAttribArray(i)
                    glVertexAttribPointer(i,4,GL_FLOAT,GL_FALSE,96,self.bondPosVBO+(i-1)*16)
                    glVertexAttribDivisor(i,1)
                self.bondPosVBO.unbind()
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.torusVBO)/3,len(self.bondPosVBO))
                for i in range(1,7):
                    glDisableVertexAttribArray(i)
                    glVertexAttribDivisor(i,0)
            else:
                for i in self.bondPosVBO:
                    self.bondShader.setUniformValue('mMatrix',QMatrix4x4(i[:16]).transposed())
                    self.bondShader.setUniformValue('s1Color',*i[16:20])
                    self.bondShader.setUniformValue('s2Color',*i[20:])
                    glDrawArrays(GL_TRIANGLES,0,len(self.torusVBO)/3)

            glDisableVertexAttribArray(0)
            self.bondShader.release()

        def drawCell(self):
                self.lineShader.bind()

                self.cellVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.cellVBO.unbind()

                self.lineShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
                self.lineShader.setUniformValue('color',QColor(0,0,0))

                glDrawArrays(GL_LINES,0,len(self.cellVBO))

                glDisableVertexAttribArray(0)
                self.lineShader.release()

        ###############################################
        # INPUT HANDLING
        ###############################################

        def mousePressEvent(self,e):
                self.setFocus()
                self.mousePos = e.pos()
                if self.mouseMode=='Modify':
                    #initiate undo
                    self.mol.init_undo()
                    #determine which atoms to modify
                    sel = self.mol.get_selection()
                    if sel:
                        self.modData=[set(i[0] for i in sel)]
                    else:
                        self.modData = [range(self.mol.get_nat())]
                    #axes in cellspace
                    mat = self.rMatrix.inverted()[0]
                    x = mat*QVector4D(1,0,0,0)
                    x = np.array([x.x(),x.y(),x.z()])
                    y = mat*QVector4D(0,-1,0,0)
                    y = np.array([y.x(),y.y(),y.z()])
                    z = mat*QVector4D(0,0,1,0)
                    z = np.array([z.x(),z.y(),z.z()])
                    self.modData+=[[x,y,z]]
                    #for rotation, center is needed
                    if e.buttons()&1:
                        vec = self.mol.get_vec()*self.mol.get_celldm()
                        def selToCoord(sel):
                            return self.mol.get_atom(sel[0])[1]+np.dot(sel[1],vec)
                        pick = self.pickAtom(e)
                        #picked atom
                        if pick:
                            self.modData+= [selToCoord(pick[1:])]
                        #com of selection
                        elif sel:
                            coords=[selToCoord(i) for i in sel]
                            self.modData+=[(np.max(coords,axis=0)+np.min(coords,axis=0))/2]
                        #com of whole cell
                        else:
                            self.modData+= [self.mol.get_center(True)]
                elif e.buttons()&2:
                    if self.mouseMode =='Camera':
                        self.rMatrix.setToIdentity()
                    elif self.mouseMode == 'Select':
                        self.mol.del_selection()
                        self.parent.updateMolStep()
                    self.updateGL()

        def mouseMoveEvent(self,e):
                delta=e.pos()-self.mousePos
                if (e.buttons() & 1):
                    if self.mouseMode == 'Camera':
                        #rotate camera
                        tmp = QMatrix4x4()
                        tmp.rotate(delta.x(),0,1,0)
                        tmp.rotate(delta.y(),1,0,0)
                        self.rMatrix = tmp*self.rMatrix
                        self.updateGL()
                    elif self.mouseMode=='Select':
                        #draw rectangle (if over threshold?)
                        pass
                    elif self.mouseMode=='Modify':
                        #rotate selected atoms around center
                        atoms=self.modData[0]
                        angle= abs(delta.x()) + abs(delta.y())
                        axis = delta.y()*self.modData[1][0]-delta.x()*self.modData[1][1]
                        shift=self.modData[2]
                        self.mol._rotate(atoms,angle,axis,shift)
                        self.parent.updateMolStep()
                        self.modMode='rotate'
                elif e.buttons()&2 and self.mouseMode=='Modify':
                    #shift selection in z-direction (cam space)
                    atoms=self.modData[0]
                    vec = self.modData[1][2]*(delta.x()+delta.y())
                    self.mol._shift(atoms,vec*0.1)
                    self.parent.updateMolStep()
                    self.modMode='shift'
                elif (e.buttons() & 4):
                    if self.mouseMode=='Camera':
                        #shift camera
                        self.xsh += delta.x()/10.
                        self.ysh -= delta.y()/10.
                        self.updateGL()
                    elif self.mouseMode=='Modify':
                        #shift selection in camera-plane
                        atoms=self.modData[0]
                        vec = self.modData[1][0]*delta.x()+self.modData[1][1]*delta.y()
                        self.mol._shift(atoms,vec*0.1)
                        self.parent.updateMolStep()
                        self.modMode='shift'
                self.mousePos = e.pos()

        def mouseReleaseEvent(self,e):
            if self.mouseMode=='Modify' and self.modMode:
                self.mol.save_undo(self.modMode)
                self.modMode=''
                self.parent.updateMolStep()
            elif self.mouseMode=='Select' and e.button()&1:
                pick = self.pickAtom(e)
                if pick:
                    self.mol.add_selection(pick[1:])
                    self.parent.updateMolStep()

        def undo(self):
            self.mol.undo()
            self.parent.updateMult()

        def pickAtom(self,coord):
            #render with selectionmode
            self.paintGL(True)

            #Wait for everything to render,configure memory alignment
            glFlush()
            glFinish()
            glPixelStorei(GL_UNPACK_ALIGNMENT,1)

            #Read pixel from GPU
            color = (GLubyte*4)(0)
            x = coord.x()
            y = self.height() - coord.y()
            glReadPixels(x,y,1,1,GL_RGBA,GL_UNSIGNED_BYTE,color)
            if color[3] == 0:
                return
            mult=self.mult[0]*self.mult[1]*self.mult[2]
            idx = color[0] + 256*color[1] + 65536*color[2]
            if idx<len(self.atomsVBO):
                realid = idx/mult
                off = idx%mult
                zoff = off%self.mult[2]
                yoff = (off/self.mult[2])%self.mult[1]
                xoff = ((off/self.mult[2])/self.mult[1])%self.mult[0]
                return [idx,realid,[xoff,yoff,zoff]]
            else:
                return None
            return idx

        def wheelEvent(self,e):
                delta = e.delta()
                #zoom with vertical wheel
                if e.orientation() & 2:
                        if delta < 0:
                                self.distance *= 1.1
                        elif delta > 0:
                                self.distance *= 0.9
                        self.updateGL()
                e.accept()

        def keyPressEvent(self,e):
                if e.key() == Qt.Key_Up:
                        tmp = QMatrix4x4()
                        tmp.rotate(-10,1,0,0)
                        self.rMatrix = tmp*self.rMatrix
                        self.updateGL()
                if e.key() == Qt.Key_Down:
                        tmp = QMatrix4x4()
                        tmp.rotate(10,1,0,0)
                        self.rMatrix = tmp*self.rMatrix
                        self.updateGL()
                if e.key() == Qt.Key_Left:
                        tmp = QMatrix4x4()
                        tmp.rotate(-10,0,1,0)
                        self.rMatrix = tmp*self.rMatrix
                        self.updateGL()
                if e.key() == Qt.Key_Right:
                        tmp = QMatrix4x4()
                        tmp.rotate(10,0,1,0)
                        self.rMatrix = tmp*self.rMatrix
                        self.updateGL()
