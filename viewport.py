#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from os.path import dirname
from math import floor
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
                self.sphereShader = QGLShaderProgram()
                self.lineShader = QGLShaderProgram()
                self.bondShader = QGLShaderProgram()
                self.selectShader = QGLShaderProgram()
                self.xsh = 0
                self.ysh = 0
                self.rMatrix = QMatrix4x4()
                self.distance = 25
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

                #set background color
                self.qglClearColor(QColor(255,255,255,0))

                #add shaders:
                self.sphereShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexSpheres.vsh')
                self.sphereShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentSpheres.fsh')
                self.bondShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexBonds.vsh')
                self.bondShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentBonds.fsh')
                self.lineShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexLines.vsh')
                self.lineShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentLines.fsh')
                self.selectShader.addShaderFromSourceFile(QGLShader.Vertex,dirname(__file__)+'/vertexSelect.vsh')
                self.selectShader.addShaderFromSourceFile(QGLShader.Fragment,dirname(__file__)+'/fragmentSelect.fsh')

                # load sphere
                sf=open(dirname(__file__)+'/sphere_model','r')
                self.sphereVBO = VBO(np.array(sf.readline().split(),'f'))
                sf.close()
                # load torus
                tf=open(dirname(__file__)+'/bond_model','r')
                self.torusVBO=VBO(np.array(tf.readline().split(),'f'))
                tf.close()

        #########################
        # Update stuff that's going to be drawn
        #########################

        def setMol(self,mol,mult):
                #prepare atoms and bonds for drawing
                if not mol:return
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
                                        tmult[i]=[x-floor(j/2) for x in range(j)]
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
                                #a = i[0]
                                #b = i[1]
                                a = atoms[i[0]][1]+i[2][0]
                                b = atoms[i[1]][1]+i[2][1]
                                #save colors
                                #TODO: remove the fucking qcolors
                                #c1 = self.pse[i[2]][2]
                                #c2 = self.pse[i[3]][2]
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

                #save atoms in VBOs
                self.atomsVBO=VBO(np.array([(at[1]+j).tolist()+[self.pse[at[0]][1],self.pse[at[0]][2].redF(),self.pse[at[0]][2].greenF(),self.pse[at[0]][2].blueF(),self.pse[at[0]][2].alphaF()] for at in atoms for j in off],'f'))

                self.updateGL()


        ##########################
        # KILL THIS SECTION WITH FIRE
        # saving the models once and for all saves startup time and grey hair
        # prepare models
        #########################
        def makeSphere(self):
                #start from octahedron:
                sphere = [QVector3D(-1.0, 0.0, 1.0),QVector3D( 1.0, 0.0, 1.0),QVector3D( 0.0, 1.0, 0.0),
                               QVector3D( 1.0, 0.0, 1.0),QVector3D( 1.0, 0.0,-1.0),QVector3D( 0.0, 1.0, 0.0),
                               QVector3D( 1.0, 0.0,-1.0),QVector3D(-1.0, 0.0,-1.0),QVector3D( 0.0, 1.0, 0.0),
                               QVector3D(-1.0, 0.0,-1.0),QVector3D(-1.0, 0.0, 1.0),QVector3D( 0.0, 1.0, 0.0),
                               QVector3D( 1.0, 0.0, 1.0),QVector3D(-1.0, 0.0, 1.0),QVector3D( 0.0,-1.0, 0.0),
                               QVector3D(-1.0, 0.0, 1.0),QVector3D(-1.0, 0.0,-1.0),QVector3D( 0.0,-1.0, 0.0),
                               QVector3D(-1.0, 0.0,-1.0),QVector3D( 1.0, 0.0,-1.0),QVector3D( 0.0,-1.0, 0.0),
                               QVector3D( 1.0, 0.0,-1.0),QVector3D( 1.0, 0.0, 1.0),QVector3D( 0.0,-1.0, 0.0)]
                #subdivide and normalize for sphere:
                for i in [0,1]:
                        sphere = self.sphere_subdiv(sphere)
                self.atom_modelspace = sphere
                #sphere2=np.array([[a.x(),a.y(),a.z()] for a in sphere],'f').flatten()
                #sf=open('sphere_model','w')
                #for i in range(len(sphere2.tolist())/3):
                #    sf.write(str(sphere2[i*3:i*3+3]).strip('[]'))
                #for i in sphere2:
                #    sf.write(str(i)+'\t')
                #sf.close()
                sf=open(dirname(__file__)+'/sphere_model','r')
                sphere2=np.array(sf.readline().split(),'f')
                sf.close()
                self.sphereVBO = VBO(sphere2)

        def sphere_subdiv(self,vertex):
                subdivide = []
                for i in range(0,len(vertex),3):
                        a = vertex[i].normalized()
                        b = vertex[i+1].normalized()
                        c = vertex[i+2].normalized()
                        d = (a+b).normalized()
                        e = (a+c).normalized()
                        f = (b+c).normalized()
                        subdivide +=[a,d,e,d,f,e,e,f,c,d,b,f]
                return subdivide

        def makeBond(self):
                #start from circle
                cylinder = [QVector3D(0.0,0.5,0.5),QVector3D(0.0,0.5,-0.5),QVector3D(0.0,-0.5,-0.5),
                            QVector3D(0.0,-0.5,0.5)]
                for i in [0,1]:
                        cylinder = self.cylinder_subdiv(cylinder)
                cylinder = self.cylinder_puzzle(cylinder)
                self.bond_modelspace = cylinder
                #write to file
                #bf=open('bond_model','w')
                #bond2=np.array([[a.x(),a.y(),a.z()] for a in cylinder],'f').flatten()
                #for i in range(len(bond2.tolist())/3):
                #    bf.write(str(bond2[i*3:i*3+3]).strip('[]'))
                #for i in bond2:
                    #bf.write(str(i)+'\t')
                #set normals
                bf=open(dirname(__file__)+'/bond_model','r')
                bond2=np.array(bf.readline().split(),'f')
                bf.close
                self.bondVBO=VBO(bond2)
                self.bo_normals_modelspace = deepcopy(cylinder)
                for i in self.bo_normals_modelspace:
                        i.setX(0)

        def cylinder_subdiv(self,vertex):
                temp = []
                for i in range(0,len(vertex)-1):
                        a = vertex[i].normalized()
                        b = vertex[i+1].normalized()
                        c = (a+b).normalized()
                        temp += [a,c]
                temp += [vertex[-1].normalized(),(vertex[-1]+vertex[0]).normalized()]
                return temp

        def cylinder_puzzle(self,vertex):
                puzzle = []
                for i in range(0,len(vertex)-1)+[-1]:
                        a = vertex[i]/3+QVector3D(-1.0,0,0)
                        b = vertex[i+1]/3+QVector3D(-1.0,0,0)
                        c = vertex[i]/3+QVector3D( 1.0,0,0)
                        d = vertex[i+1]/3+QVector3D( 1.0,0,0)
                        e = vertex[i]/3+QVector3D( 1.0,0,0)
                        f = vertex[i+1]/3+QVector3D(-1.0,0,0)
                        puzzle += [a,c,b,d,f,e]
                return puzzle

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

                #don't do anything if there's no molecule
                #if not hasattr(self,'mol'):return

                #initialize transformation matrices:
                self.mMatrix = QMatrix4x4()
                self.vMatrix = QMatrix4x4()

                #perspective
                self.vMatrix.lookAt(QVector3D(0,0,self.distance),QVector3D(0,0,0),QVector3D(0,1,0))

                #translate the view point:
                self.vMatrix.translate(self.xsh,self.ysh,0)

                #TODO: orthogonal zooming needs fix
                #check for projection:
                if self.parent.view:
                        self.proj = self.pMatrix
                else:
                        self.proj = self.oMatrix
                        #scale based on distance for zoom effect
                        self.vMatrix.scale(10./self.distance)
                #TODO: decrease quality with increasing number of atoms
                #rendering:
                if not hasattr(self,'atomsVBO'): return
                if select:
                        self.drawAtomsSelect()
                else:
                        self.drawAtoms()
                        #self.drawAtomsSelect()
                        if self.parent.bondShow:
                                self.drawBonds()
                        else:
                                pass
                        if self.parent.cell:
                                self.drawCell()

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

        def drawAtomsSelect(self):
                #bind shaders:
                self.selectShader.bind()

                #send vertices
                self.sphereVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.sphereVBO.unbind()

                # transformation matrices, model matrix needs not perform anything here
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

        def drawBonds(self):
                self.bondShader.bind()

                #send vertices
                self.torusVBO.bind()
                glEnableVertexAttribArray(0)
                glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
                self.torusVBO.unbind()

                #send positions and rotations via mMatrix
                self.bondPosVBO.bind()
                for i in range(1,7):
                    glEnableVertexAttribArray(i)
                    glVertexAttribPointer(i,4,GL_FLOAT,GL_FALSE,96,self.bondPosVBO+(i-1)*16)
                    glVertexAttribDivisor(i,1)
                self.bondPosVBO.unbind()

                #load model Matrix with coordinates relative to offset
                self.mMatrix.setToIdentity()
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


        def drawSelection(self):
                #bind shaders:
                self.sphereShader.bind()

                #send vertices
                self.sphereShader.setAttributeArray('vertex_modelspace',self.atom_modelspace)
                self.sphereShader.enableAttributeArray('vertex_modelspace')
                #send normals for lighting
                #equal coordinates for sphere
                self.sphereShader.setAttributeArray('normals_modelspace',self.atom_modelspace)
                self.sphereShader.enableAttributeArray('normals_modelspace')

                #render loop
                for i in self.parent.selection:
                        #load model matrix with coordinates
                        self.mMatrix.setToIdentity()
                        #move atoms to coord and recognize offset
                        atom = self.parent.atoms[i]
                        self.mMatrix.translate(atom[0][0],atom[0][1],atom[0][2])
                        #scale atoms depending on species
                        self.mMatrix.scale(atom[2]*1.3)
                        #bind transformation matrices
                        self.sphereShader.setUniformValue('mvpMatrix',self.proj*self.vMatrix*self.mMatrix)
                        self.sphereShader.setUniformValue('mvMatrix',self.vMatrix*self.mMatrix)
                        #color vertices
                        self.sphereShader.setUniformValue('MaterialDiffuseColor',QColor(120,120,150,200))
                        #glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
                        #draw
                        glDrawArrays(GL_TRIANGLES,0,len(self.atom_modelspace))
                #reset
                self.sphereShader.disableAttributeArray('vertex_modelspace')
                self.sphereShader.disableAttributeArray('normals_modelspace')
                self.sphereShader.release()

        def drawCenter(self):
                self.lineShader.bind()
                self.mMatrix.setToIdentity()
                self.lineShader.setUniformValue('mvpMatrix',self.proj*self.vMatrix*self.mMatrix)
                self.lineShader.setUniformValue('color',QColor(0,0,0))
                self.lineShader.setAttributeArray('vertex_modelspace',self.atom_modelspace)
                self.lineShader.enableAttributeArray('vertex_modelspace')
                glDrawArrays(GL_TRIANGLES,0,len(self.atom_modelspace))
                self.lineShader.disableAttributeArray('vertex_modelspace')
                self.lineShader.release()

        ###############################################
        # INPUT HANDLING
        ###############################################

        def mousePressEvent(self,e):
                if self.parent.mouseMode!=True:return
                self.setFocus()
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
                if self.parent.mouseMode != True:return
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
                if self.parent.mouseMode != False:return
                #here be render function
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
                if color[3]==0:
                        self.parent.selection=[]
                else:
                        id = color[0] + 256*color[1] + 65536*color[2]
                        ctrl = e.modifiers() == Qt.ControlModifier
                        prev = id in self.parent.selection
                        if ctrl and prev:
                                del self.parent.selection[self.parent.selection.index(id)]
                        elif ctrl and not prev:
                                self.parent.selection.append(id)
                        elif prev and not ctrl and len(self.parent.selection)>1:
                                self.parent.selection = [id]
                        elif prev and not ctrl and len(self.parent.selection)<2:
                                self.parent.selection = []
                        elif not prev and not ctrl:
                                self.parent.selection = [id]
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
