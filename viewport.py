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
                form.setVersion(3,3)
                form.setProfile(QGLFormat.CoreProfile)
                super(ViewPort,self).__init__(form)
                self.parent = parent
                self.showBonds = True
                self.showCell = True
                self.showPlane = False
                self.showSurf = False
                self.planeVal = 0
                self.perspective = True
                self.mouseSelect = False
                self.AA = True
                self.xsh = 0
                self.ysh = 0
                self.sel = []
                self.rMatrix = QMatrix4x4()
                self.distance = 25
                self.initActions()
                # radii, color
                self.pse={'H' : [1.20,0.38,QColor(191,191,191,255)],
                          'He': [1.40,0.32,QColor(216,255,255,255)],
                          'Li': [1.82,1.34,QColor(204,127,255,255)],
                          'Be': [1.53,0.90,QColor(193,255,  0,255)],
                          'B' : [1.92,0.82,QColor(255,181,181,255)],
                          'C' : [1.70,0.77,QColor(102,102,102,255)],
                          'N' : [1.55,0.75,QColor( 12, 12,255,255)],
                          'O' : [1.52,0.73,QColor(255, 12, 12,255)],
                          'F' : [1.47,0.71,QColor(127,178,255,255)],
                          'Ne': [1.54,0.69,QColor(178,226,244,255)],
                          'Na': [2.27,1.54,QColor(170, 91,242,255)],
                          'Mg': [1.73,1.30,QColor(137,255,  0,255)],
                          'Al': [1.84,1.18,QColor(191,165,165,255)],
                          'Si': [2.10,1.11,QColor(127,153,153,255)],
                          'P' : [1.80,1.06,QColor(255,127,  0,255)],
                          'S' : [1.80,1.02,QColor(178,178,  0,255)],
                          'Cl': [1.75,0.99,QColor( 30,239, 30,255)],
                          'Ar': [1.88,0.97,QColor(127,209,226,255)],
                          'K' : [2.75,1.96,QColor(142, 63,211,255)],
                          'Ca': [2.31,1.74,QColor( 61,255,  0,255)],
                          'Sc': [2.11,1.44,QColor(229,229,229,255)],
                          'Ti': [1.70,1.36,QColor(191,193,198,255)],
                          'V' : [1.70,1.25,QColor(165,165,170,255)],
                          'Cr': [1.70,1.27,QColor(137,153,198,255)],
                          'Mn': [1.70,1.39,QColor(155,122,198,255)],
                          'Fe': [1.70,1.25,QColor(224,102, 51,255)],
                          'Co': [1.70,1.26,QColor(239,142,160,255)],
                          'Ni': [1.63,1.21,QColor( 79,209, 79,255)],
                          'Cu': [1.40,1.38,QColor(198,127, 51,255)],
                          'Zn': [1.39,1.31,QColor(124,127,175,255)],
                          'Ga': [1.87,1.26,QColor(193,142,142,255)],
                          'Ge': [2.11,1.22,QColor(102,142,142,255)],
                          'As': [1.85,1.19,QColor(188,127,226,255)],
                          'Se': [1.90,1.16,QColor(255,160,  0,255)],
                          'Br': [1.85,1.14,QColor(165, 40, 40,255)],
                          'Kr': [2.02,1.10,QColor( 91,183,209,255)],
                          'Rb': [3.03,2.11,QColor(112, 45,175,255)],
                          'Sr': [2.49,1.92,QColor(  0,255,  0,255)],
                          'Y' : [1.70,1.62,QColor(147,255,255,255)],
                          'Zr': [1.70,1.48,QColor(147,224,224,255)],
                          'Nb': [1.70,1.37,QColor(114,193,201,255)],
                          'Mo': [1.70,1.45,QColor( 84,181,181,255)],
                          'Tc': [1.70,1.56,QColor( 58,158,158,255)],
                          'Ru': [1.70,1.26,QColor( 35,142,142,255)],
                          'Rh': [1.70,1.35,QColor( 10,124,140,255)],
                          'Pd': [1.63,1.31,QColor(  0,104,132,255)],
                          'Ag': [1.72,1.53,QColor(224,224,255,255)],
                          'Cd': [1.58,1.48,QColor(255,216,142,255)],
                          'In': [1.93,1.44,QColor(165,117,114,255)],
                          'Sn': [2.17,1.41,QColor(102,127,127,255)],
                          'Sb': [2.06,1.38,QColor(158,99 ,181,255)],
                          'Te': [2.06,1.35,QColor(211,122,  0,255)],
                          'I' : [1.98,1.33,QColor(147,  0,147,255)],
                          'Xe': [2.16,1.30,QColor( 66,158,175,255)],
                          'Cs': [3.43,2.25,QColor( 86, 22,142,255)],
                          'Ba': [2.68,1.98,QColor(  0,201,  0,255)],
                          'La': [1.70,1.69,QColor(112,211,255,255)],
                          'Ce': [1.70,0.77,QColor(255,255,198,255)],
                          'Pr': [1.70,0.77,QColor(216,255,198,255)],
                          'Nd': [1.70,0.77,QColor(198,255,198,255)],
                          'Pm': [1.70,0.77,QColor(163,255,198,255)],
                          'Sm': [1.70,0.77,QColor(142,255,198,255)],
                          'Eu': [1.70,0.77,QColor( 96,255,198,255)],
                          'Gd': [1.70,0.77,QColor( 68,255,198,255)],
                          'Tb': [1.70,0.77,QColor( 48,255,198,255)],
                          'Dy': [1.70,0.77,QColor( 30,255,198,255)],
                          'Ho': [1.70,0.77,QColor(  0,255,155,255)],
                          'Er': [1.70,0.77,QColor(  0,229,117,255)],
                          'Tm': [1.70,0.77,QColor(  0,211, 81,255)],
                          'Yb': [1.70,0.77,QColor(  0,191, 56,255)],
                          'Lu': [1.70,1.60,QColor(  0,170, 35,255)],
                          'Hf': [1.70,1.50,QColor( 76,193,255,255)],
                          'Ta': [1.70,1.38,QColor( 76,165,255,255)],
                          'W' : [1.70,1.46,QColor( 33,147,214,255)],
                          'Re': [1.70,1.59,QColor( 38,124,170,255)],
                          'Os': [1.70,1.28,QColor( 38,102,150,255)],
                          'Ir': [1.70,1.37,QColor( 22, 84,135,255)],
                          'Pt': [1.75,1.28,QColor(244,237,209,255)],
                          'Au': [1.66,1.44,QColor(204,209, 30,255)],
                          'Hg': [1.55,1.49,QColor(181,181,193,255)],
                          'Tl': [1.96,1.48,QColor(165, 84, 76,255)],
                          'Pb': [2.02,1.47,QColor( 86, 89, 96,255)],
                          'Bi': [2.07,1.46,QColor(158, 79,181,255)],
                          'Po': [1.97,0.77,QColor(170, 91,  0,255)],
                          'At': [2.02,0.77,QColor(117, 79, 68,255)],
                          'Rn': [2.20,1.45,QColor( 66,130,150,255)],
                          'Fr': [3.48,0.77,QColor( 66,  0,102,255)],
                          'Ra': [2.83,0.77,QColor(0  ,124,  0,255)],
                          'Ac': [1.70,0.77,QColor(112,170,249,255)],
                          'Th': [1.70,0.77,QColor(0  ,186,255,255)],
                          'Pa': [1.70,0.77,QColor(0  ,160,255,255)],
                          'U' : [1.86,0.77,QColor(0  ,142,255,255)],
                          'Np': [1.70,0.77,QColor(0  ,127,255,255)],
                          'Pu': [1.70,0.77,QColor(0  ,107,255,255)],
                          'Am': [1.70,0.77,QColor(84 , 91,242,255)],
                          'Cm': [1.70,0.77,QColor(119, 91,226,255)],
                          'Bk': [1.70,0.77,QColor(137, 79,226,255)],
                          'Cf': [1.70,0.77,QColor(160, 53,211,255)],
                          'Es': [1.70,0.77,QColor(178, 30,211,255)],
                          'Fm': [1.70,0.77,QColor(178, 30,186,255)],
                          'Md': [1.70,0.77,QColor(178, 12,165,255)],
                          'No': [1.70,0.77,QColor(188, 12,135,255)],
                          'Lr': [1.70,0.77,QColor(198,  0,102,255)],
                          'Rf': [1.70,0.77,QColor(204,  0, 89,255)],
                          'Db': [1.70,0.77,QColor(209,  0, 79,255)],
                          'Sg': [1.70,0.77,QColor(216,  0, 68,255)],
                          'Bh': [1.70,0.77,QColor(224,  0, 56,255)],
                          'Hs': [1.70,0.77,QColor(229,  0, 45,255)],
                          'Mt': [1.70,0.77,QColor(234,  0, 38,255)],
                          'Ds': [1.70,0.77,QColor(237,  0, 35,255)],
                          'Rg': [1.70,0.77,QColor(239,  0, 33,255)],
                          'Cn': [1.70,0.77,QColor(242,  0, 30,255)],
                          'Uut':[1.70,0.77,QColor(244,  0, 28,255)],
                          'Fl': [1.70,0.77,QColor(247,  0, 25,255)],
                          'Uup':[1.70,0.77,QColor(249,  0, 22,255)],
                          'Lv': [1.70,0.77,QColor(252,  0, 20,255)],
                          'Uus':[1.70,0.77,QColor(252,  0, 17,255)],
                          'Uuo':[1.70,0.77,QColor(252,  0, 15,255)]}


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

                #add shaders:
                self.sphereShader = QGLShaderProgram()
                self.sphereShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexSpheres.vsh')
                self.sphereShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentSpheres.fsh')

                self.bondShader = QGLShaderProgram()
                self.bondShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexBonds.vsh')
                self.bondShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentBonds.fsh')

                self.lineShader = QGLShaderProgram()
                self.lineShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexLines.vsh')
                self.lineShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentLines.fsh')

                self.selectShader = QGLShaderProgram()
                self.selectShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexSelect.vsh')
                self.selectShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentSelect.fsh')

                self.planeShader = QGLShaderProgram()
                self.planeShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexPlane.vsh')
                self.planeShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentPlane.fsh')
                self.surfShader = QGLShaderProgram()
                self.surfShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexSurf.vsh')
                self.surfShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentSurf.fsh')

                # load sphere
                sf=open(dirname(__file__)+'/sphere_model','r')
                self.sphereVBO = VBO(np.array(sf.readline().split(),'f'))
                sf.close()
                # load torus
                tf=open(dirname(__file__)+'/bond_model','r')
                self.torusVBO=VBO(np.array(tf.readline().split(),'f'))
                tf.close()

        ##########################################
        # Modify state of Visualization
        ##########################################
        def initActions(self):
                changeProj = QAction('Change &Projection',self)
                changeProj.setShortcut('p')
                changeProj.triggered.connect(self.changeState)
                toggleCell = QAction('Toggle &Cell',self)
                toggleCell.setShortcut('c')
                toggleCell.triggered.connect(self.changeState)
                mouseRotate = QAction('&Rotate',self)
                mouseRotate.setShortcut('r')
                mouseRotate.triggered.connect(self.changeState)
                mouseSelect = QAction('&Select',self)
                mouseSelect.setShortcut('s')
                mouseSelect.triggered.connect(self.changeState)
                antiAlias = QAction('Toggle &Antialiasing',self)
                antiAlias.setShortcut('a')
                antiAlias.triggered.connect(self.changeState)
                toggleBonds = QAction('Toggle &Bonds',self)
                toggleBonds.setShortcut('b')
                toggleBonds.triggered.connect(self.changeState)
                vMenu = self.parent.parent.menuBar().addMenu('&View')
                vMenu.addAction(changeProj)
                vMenu.addAction(toggleCell)
                vMenu.addAction(toggleBonds)
                vMenu.addAction(antiAlias)
                vMenu.addSeparator()
                vMenu.addAction(mouseRotate)
                vMenu.addAction(mouseSelect)

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
                elif reason=='&Rotate':
                    self.mouseSelect=False
                elif reason=='&Select':
                    self.mouseSelect=True
                elif reason=='Toggle &Antialiasing':
                    if self.AA:
                            glDisable(GL_MULTISAMPLE)
                    elif not self.AA:
                            glEnable(GL_MULTISAMPLE)
                    self.AA = not self.AA
                    self.updateGL()

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

                #clear old selection if present:
                if hasattr(self,'selVBO'):
                    del self.selVBO

                #save for interaction
                self.mol=mol
                self.mult=mult
                #local variables for convenience
                atoms = mol.get_all_atoms()
                vec = mol.get_vec()*mol.get_celldm()
                center = mol.get_center()
                bonds = mol.get_bonds()

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
                tempbonds=[]
                #binary representation of enabled pbc directions (b'zyx')
                mult=np.sign(mult[0]-1)+np.sign(mult[1]-1)*2+np.sign(mult[2]-1)*4
                for j in [0,1,2,3,4,5,6,7]:
                        #check if pbc part is necessary
                        if mult&j!=j: continue
                        for i in bonds[j]:
                                #get positions of atoms
                                a = atoms[i[0]][1]+i[2][0]
                                b = atoms[i[1]][1]+i[2][1]
                                #save colors
                                c1 = self.pse[atoms[i[0]][0]][2]
                                c2 = self.pse[atoms[i[1]][0]][2]
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
                                        tempbonds+=[ic*ax[0]*ax[0]+c,ic*ax[0]*ax[1]-s*ax[2],ic*ax[0]*ax[2]+s*ax[1],0.,
                                                    ic*ax[0]*ax[1]+s*ax[2],ic*ax[1]*ax[1]+c,ic*ax[1]*ax[2]-s*ax[0],0.,
                                                    ic*ax[0]*ax[2]-s*ax[1],ic*ax[1]*ax[2]+s*ax[0],ic*ax[2]*ax[2]+c,0.,
                                                    pos[0]+k[0],pos[1]+k[1],pos[2]+k[2],1.,
                                                    c1.redF(),c1.greenF(),c1.blueF(),c1.alphaF(),
                                                    c2.redF(),c2.greenF(),c2.blueF(),c2.alphaF()]
                self.bondPosVBO=VBO(np.array(tempbonds,'f'))

                #save atoms in VBOs
                self.atomsVBO=VBO(np.array([(at[1]+j).tolist()+[self.pse[at[0]][1],self.pse[at[0]][2].redF(),self.pse[at[0]][2].greenF(),self.pse[at[0]][2].blueF(),self.pse[at[0]][2].alphaF()] for at in atoms for j in off],'f'))
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
                self.offVBO=VBO(np.array(off,'f'))

                self.updateGL()

        ################################################
        # CREATE AND MANAGE PLANES AND SURFACES
        ################################################
        def toggleSurf(self):
            self.showSurf = not self.showSurf
            self.updateGL()

        def setSurf(self,sval):
            vertices,nv = make_iso_surf(self.mol.get_vol(),sval)
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
                #bind shaders:
                self.sphereShader.bind()

                #send vertices
                self.sphereVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.sphereVBO.unbind()

                # transformation matrices, model matrix needs not perform anything here
                self.sphereShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
                self.sphereShader.setUniformValue('rMatrix',self.rMatrix)

                #interleaved VBO
                self.atomsVBO.bind()
                #position vector
                glEnableVertexAttribArray(1)
                glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
                #scale factors
                glEnableVertexAttribArray(2)
                glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.atomsVBO+12)
                #color vectors
                glEnableVertexAttribArray(3)
                glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,32,self.atomsVBO+16)
                glVertexAttribDivisor(1,1)
                glVertexAttribDivisor(2,1)
                glVertexAttribDivisor(3,1)
                self.atomsVBO.unbind()
                # draw instances
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO),len(self.atomsVBO))

                #reset
                glDisableVertexAttribArray(0)
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glDisableVertexAttribArray(3)
                glVertexAttribDivisor(1,0)
                glVertexAttribDivisor(2,0)
                glVertexAttribDivisor(3,0)
                self.sphereShader.release()

        def drawSurf(self):
            self.surfShader.bind()

            #send vertices
            self.surfVBO.bind()
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
            self.surfVBO.unbind()

            self.surfShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
            self.surfShader.setUniformValue('cellVec',QMatrix3x3((self.mol.get_vec()*self.mol.get_celldm()).flatten()))

            self.offVBO.bind()
            glEnableVertexAttribArray(1)
            glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,None)
            glVertexAttribDivisor(1,1)
            self.offVBO.unbind()

            #glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
            glDisable(GL_CULL_FACE)
            glDrawArraysInstanced(GL_TRIANGLES,0,len(self.surfVBO),len(self.offVBO))
            #glDrawArraysInstanced(GL_POINTS,0,np.prod(shape),len(self.offVBO))
            glEnable(GL_CULL_FACE)
            #glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)

            glDisableVertexAttribArray(0)
            glDisableVertexAttribArray(1)
            glVertexAttribDivisor(1,0)
            self.surfShader.release()

        def drawPlane(self):
                self.planeShader.bind()

                #send plane vertices
                self.planeVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,20,None)
                glEnableVertexAttribArray(1)
                glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,20,self.planeVBO+12)
                self.planeVBO.unbind()

                #offset for instances
                self.offVBO.bind()
                glEnableVertexAttribArray(2)
                glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,0,None)
                glVertexAttribDivisor(2,1)
                self.offVBO.unbind()

                #transformation matrices
                self.planeShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)

                #send texture
                ID = glGenTextures(1)
                glActiveTexture(GL_TEXTURE0)
                glBindTexture(GL_TEXTURE_2D,ID)
                glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
                glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
                glTexImage2D(GL_TEXTURE_2D,0,GL_RED,self.planeTex.shape[1],self.planeTex.shape[0],0,GL_RED,GL_FLOAT,self.planeTex)
                self.planeShader.setUniformValue('texSampler',0)

                #render with disabled culling
                glDisable(GL_CULL_FACE)
                glDrawArraysInstanced(GL_TRIANGLE_STRIP,0,len(self.planeVBO),len(self.offVBO))
                glEnable(GL_CULL_FACE)

                #reset
                glDeleteTextures(1)
                glDisableVertexAttribArray(0)
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glVertexAttribDivisor(2,0)
                self.planeShader.release()

        def drawAtomsSelect(self):
                #bind shaders:
                self.selectShader.bind()

                #send vertices
                self.sphereVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.sphereVBO.unbind()

                # transformation matrices
                self.selectShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)

                #interleaved VBO
                self.atomsVBO.bind()
                #position vector
                glEnableVertexAttribArray(1)
                glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
                #scale factors
                glEnableVertexAttribArray(2)
                glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.atomsVBO+12)
                #colors are based on gl_InstanceID in vertex shader
                glVertexAttribDivisor(1,1)
                glVertexAttribDivisor(2,1)
                self.atomsVBO.unbind()
                # draw instances
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO),len(self.atomsVBO))

                #reset
                glDisableVertexAttribArray(0)
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glVertexAttribDivisor(1,0)
                glVertexAttribDivisor(2,0)
                self.selectShader.release()

        def drawSelection(self):
                #bind shaders:
                self.sphereShader.bind()

                #send vertices
                self.sphereVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.sphereVBO.unbind()

                # transformation matrices, model matrix needs not perform anything here
                self.sphereShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
                self.sphereShader.setUniformValue('rMatrix',self.rMatrix)

                #interleaved VBO
                self.selVBO.bind()
                #position vector
                glEnableVertexAttribArray(1)
                glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
                #scale factors
                glEnableVertexAttribArray(2)
                glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.selVBO+12)
                #color vectors
                glEnableVertexAttribArray(3)
                glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,32,self.selVBO+16)
                glVertexAttribDivisor(1,1)
                glVertexAttribDivisor(2,1)
                glVertexAttribDivisor(3,1)
                self.selVBO.unbind()
                # draw instances
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO),len(self.selVBO))

                #reset
                glDisableVertexAttribArray(0)
                glDisableVertexAttribArray(1)
                glDisableVertexAttribArray(2)
                glDisableVertexAttribArray(3)
                glVertexAttribDivisor(1,0)
                glVertexAttribDivisor(2,0)
                glVertexAttribDivisor(3,0)
                self.sphereShader.release()

        def drawBonds(self):
                self.bondShader.bind()

                #send vertices
                self.torusVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.torusVBO.unbind()

                #send positions and rotations via interleaved VBO
                self.bondPosVBO.bind()
                for i in range(1,7):
                    glEnableVertexAttribArray(i)
                    glVertexAttribPointer(i,4,GL_FLOAT,GL_FALSE,96,self.bondPosVBO+(i-1)*16)
                    glVertexAttribDivisor(i,1)
                self.bondPosVBO.unbind()

                # bind transformation matrices
                self.bondShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
                self.bondShader.setUniformValue('rMatrix',self.rMatrix)
                #color vertices
                self.bondShader.setUniformValue('Side1DiffuseColor',self.pse['C'][2])
                self.bondShader.setUniformValue('Side2DiffuseColor',self.pse['C'][2])
                #draw
                glDrawArraysInstanced(GL_TRIANGLES,0,len(self.torusVBO.data),len(self.bondPosVBO)/24)

                #reset
                for i in range(1,7):
                    glDisableVertexAttribArray(i)
                    glVertexAttribDivisor(i,0)
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
                if self.mouseSelect:return
                if (e.buttons() & 1):
                        self.oldX = self.newX = e.x()
                        self.oldY = self.newY = e.y()
                elif (e.buttons() & 2):
                        self.rMatrix.setToIdentity()
                self.updateGL()
                #store initial position
                self.mousePos = e.pos()
                #stop event from getting pushed to parent
                e.accept()

        def mouseMoveEvent(self,e):
                if self.mouseSelect:return
                #calculate deviation from initial position
                deltaX = e.x() - self.mousePos.x()
                deltaY = e.y() - self.mousePos.y()

                #left click: rotate molecule
                if (e.buttons() & 1):
                        #get new position
                        self.newX = e.x()
                        self.newY = e.y()
                        #apply rotation and store it (important!)
                        tmp = QMatrix4x4()
                        deltaX = self.newX-self.oldX
                        deltaY = self.newY - self.oldY
                        tmp.rotate(deltaX,0,1,0)
                        tmp.rotate(deltaY,1,0,0)
                        self.rMatrix = tmp*self.rMatrix
                        #save as old positions for next step
                        self.oldX = e.x()
                        self.oldY = e.y()
                        self.updateGL()
                elif (e.buttons() & 4):
                        self.xsh += deltaX/10.
                        self.ysh -= deltaY/10.
                        self.updateGL()

                self.mousePos = e.pos()
                e.accept()

        def mouseReleaseEvent(self,e):
                if not self.mouseSelect:return
                if(e.button() & 2):
                    self.sel=[]
                    if hasattr(self,'selVBO'):
                        del self.selVBO
                elif(e.button()&4) and self.sel:
                    self.sel.pop()
                    if hasattr(self,'selVBO'):
                        if len(self.sel)==0:
                            del self.selVBO
                        else:
                            self.selVBO=VBO(np.array([self.atomsVBO.data[a][:3].tolist()+
                                [self.atomsVBO.data[a][3]*1.5,
                                    0.4,0.4,0.5,0.5] for a in [b[0] for b in self.sel]],'f'))
                elif(e.button()&1) and len(self.sel)<4:
                    #render with selectionmode
                    self.paintGL(True)

                    #Wait for everything to render,configure memory alignment
                    glFlush()
                    glFinish()
                    glPixelStorei(GL_UNPACK_ALIGNMENT,1)

                    #Read pixel from GPU
                    color = (GLubyte*4)(0)
                    x = e.x()
                    y = self.height() - e.y()
                    glReadPixels(x,y,1,1,GL_RGBA,GL_UNSIGNED_BYTE,color)
                    if color[3] == 0:
                        return
                    mult=self.mult[0]*self.mult[1]*self.mult[2]
                    id = color[0] + 256*color[1] + 65536*color[2]
                    if id<len(self.atomsVBO.data):
                        if id in [a[0] for a in self.sel]:
                            return
                        else:
                            at = self.mol.get_atom(id/mult)
                            self.sel.append([id,id/mult,at[0],
                                    self.atomsVBO.data[id][0:3]])
                            self.selVBO=VBO(np.array([self.atomsVBO.data[a][:3].tolist()+
                                [self.atomsVBO.data[a][3]*1.5,
                                    0.4,0.4,0.5,0.5] for a in [b[0] for b in self.sel]],'f'))
                    else:
                        return
                self.parent.edit.setSel(self.sel)
                self.updateGL()

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
