# -*- coding: utf-8 -*-

from os.path import dirname
from copy import deepcopy
import numpy as np

from PyQt4.QtCore import Qt,QRectF,QTimer
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
from OpenGL.GL import *
from OpenGL.GL.shaders import *
from OpenGL.arrays.vbo import *

from .gui_c import makeIsoSurf
from ..settings import config

class VisualWidget(QFrame):

    def __init__(self,parent):
        super(VisualWidget,self).__init__()
        self.parent = parent
        self.viewport = ViewPort(parent)
        self.updatedisable = False
        mouseLayout = self.initMouse()
        control1 = self.initMod()
        control2 = self.initControl()
        vbox = QVBoxLayout()
        vbox.addLayout(mouseLayout)
        vbox.addWidget(self.viewport)
        vbox.addLayout(control1)
        vbox.addLayout(control2)
        vbox.setStretchFactor(self.viewport,1)
        vbox.setContentsMargins(0,0,0,0)
        self.setLayout(vbox)
        self.setFrameStyle(38)

    def initMouse(self):
        self.mouseGroup = QButtonGroup()
        mouseLayout = QHBoxLayout()
        camBut = QPushButton("Camera")
        camBut.setShortcut("r")
        camBut.setToolTip("Move camera\n\nLMB: Rotate camera\nMMB: Drag camera\nRMB: Align camera to z-axis")
        selBut = QPushButton("Select")
        selBut.setShortcut("s")
        selBut.setToolTip("Select atoms\n\nLMB: Select atoms\nRMB: Clear selection")
        modBut = QPushButton("Modify")
        modBut.setShortcut("m")
        modBut.setToolTip("Modify geometry\n\nLMB: Rotate atoms (around Center of Mass or selected Atom)\nMMB: Move atoms in xy-plane (camera)\nRMB: Move atoms along z-axis (camera)")
        for i in [camBut,selBut,modBut]:
            mouseLayout.addWidget(i)
            self.mouseGroup.addButton(i)
            i.setCheckable(True)
        self.mouseGroup.buttonClicked.connect(self.viewport.setMouseMode)
        self.mouseGroup.setExclusive(True)
        camBut.setChecked(True)
        return mouseLayout

    def initMod(self):
        row1 = QHBoxLayout()
        #Camera-Alignment buttons
        row1.addWidget(QLabel("Align camera:"))
        for i,w in enumerate(["+x","+y","+z","-x","-y","-z"]):
            button = QPushButton(w)
            button.clicked.connect(self.viewport.alignView)
            button.setFixedWidth(30)
            row1.addWidget(button)
        row1.addStretch()
        #Cell multiplication
        self.mult = []
        row1.addWidget(QLabel("Cell multiply:"))
        for i,s in enumerate(["x:","y:","z:"]):
            row1.addWidget(QLabel(s))
            self.mult.append(QSpinBox())
            self.mult[-1].setMinimum(1)
            self.mult[-1].valueChanged.connect(self.updateView)
            row1.addWidget(self.mult[-1])
        return row1

    def initControl(self):
        ##Choose step, animate
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Select step: "))
        self.buttons = []
        for i in range(5):
            self.buttons.append(QPushButton())
            self.buttons[-1].clicked.connect(self.changeStep)
        self.buttons[0].setIcon(self.style().standardIcon(64))
        self.buttons[0].setFixedWidth(30)
        self.buttons[1].setAutoRepeat(True)
        self.buttons[1].setIcon(self.style().standardIcon(66))
        self.buttons[1].setFixedWidth(30)
        self.buttons[2].setIcon(self.style().standardIcon(60))
        self.buttons[2].setFixedWidth(30)
        self.buttons[3].setIcon(self.style().standardIcon(65))
        self.buttons[3].setAutoRepeat(True)
        self.buttons[3].setFixedWidth(30)
        self.buttons[4].setIcon(self.style().standardIcon(63))
        self.buttons[4].setFixedWidth(30)
        self.animTimer = QTimer()
        self.animTimer.setInterval(50)
        self.animTimer.timeout.connect(self.buttons[3].click)
        for i,w in enumerate(self.buttons):
            row2.addWidget(w)
        self.slideStep = QSlider()
        self.slideStep.setOrientation(1)
        self.slideStep.setMinimum(1)
        self.slideStep.setMaximum(1)
        self.slideStep.setTickPosition(self.slideStep.TicksBelow)
        self.slideStep.setSingleStep(1)
        self.slideStep.valueChanged.connect(self.changeStep)
        row2.addWidget(self.slideStep)
        self.stepRange = QIntValidator()
        self.stepRange.setBottom(1)
        self.textStep = QLineEdit("1")
        self.textStep.editingFinished.connect(self.changeStep)
        self.textStep.setValidator(self.stepRange)
        self.textStep.setFixedWidth(30)
        row2.addWidget(self.textStep)
        row2.addWidget(QLabel("/"))
        self.maxStep = QLabel("1")
        row2.addWidget(self.maxStep)
        return row2

    def setMol(self,mol):
        self.mol = mol
        step = mol.curStep+1
        length = len(mol)
        if length == 1:
            self.textStep.setDisabled(True)
            self.slideStep.setDisabled(True)
        else:
            self.textStep.setEnabled(True)
            self.slideStep.setEnabled(True)
        if step == 1:
            self.buttons[0].setDisabled(True)
            self.buttons[1].setDisabled(True)
        else:
            self.buttons[0].setEnabled(True)
            self.buttons[1].setEnabled(True)
        if step == length:
            self.animTimer.stop()
            self.buttons[3].setDisabled(True)
            self.buttons[4].setDisabled(True)
        else:
            self.buttons[3].setEnabled(True)
            self.buttons[4].setEnabled(True)
        if length and step<length:
            self.buttons[2].setEnabled(True)
        else:
            self.buttons[2].setDisabled(True)
        self.stepRange.setTop(length)
        self.maxStep.setText(str(length))
        self.textStep.setText(str(step))
        self.slideStep.blockSignals(True)
        self.slideStep.setMaximum(length)
        self.slideStep.setValue(step)
        self.slideStep.blockSignals(False)
        self.updateView()

    def updateView(self):
        self.viewport.setMol(self.mol,[i.value() for i in self.mult])

    def changeStep(self):
        sender = self.sender()
        if sender is self.slideStep:
            step = self.slideStep.value()-1
            self.textStep.setText(str(step+1))
        elif sender is self.textStep:
            step = int(self.textStep.text())-1
            self.slideStep.blockSignals(True)
            self.slideStep.setValue(step+1)
            self.slideStep.blockSignals(False)
        elif sender in self.buttons:
            index = self.buttons.index(sender)
            if index == 0:
                step = 0
            elif index == 1:
                step = self.slideStep.value()-2
            elif index == 2:
                if self.animTimer.isActive():
                        self.animTimer.stop()
                else:
                        self.animTimer.start()
                return
            elif index == 3:
                step = self.slideStep.value()
            elif index == 4:
                step = len(self.mol)-1
        self.mol.changeStep(step)
        self.parent.updateMol()

class ViewPort(QGLWidget):

    def __init__(self,parent):
        #init with antialiasing, transparency and OGLv3.3 core profile
        form=QGLFormat(QGL.SampleBuffers|QGL.AlphaChannel)
        form.setProfile(QGLFormat.CoreProfile)
        super(ViewPort,self).__init__(form)
        self.parent = parent
        self.mouseMode = "Camera"
        self.mousePos = None
        self.rectPos = None
        self.xsh = 0
        self.ysh = 0
        self.rMatrix = QMatrix4x4()
        self.distance = 25
        self.copyBuf=[]

    ###############################################
    # INPUT HANDLING
    ###############################################

    def keyPressEvent(self,e):
        if e.key() == Qt.Key_Up:
                tmp = QMatrix4x4()
                tmp.rotate(-10,1,0,0)
                self.rMatrix = tmp*self.rMatrix
                self.update()
        if e.key() == Qt.Key_Down:
                tmp = QMatrix4x4()
                tmp.rotate(10,1,0,0)
                self.rMatrix = tmp*self.rMatrix
                self.update()
        if e.key() == Qt.Key_Left:
                tmp = QMatrix4x4()
                tmp.rotate(-10,0,1,0)
                self.rMatrix = tmp*self.rMatrix
                self.update()
        if e.key() == Qt.Key_Right:
                tmp = QMatrix4x4()
                tmp.rotate(10,0,1,0)
                self.rMatrix = tmp*self.rMatrix
                self.update()

    def setMouseMode(self,but):
        t=but.text()
        self.mouseMode = t
        if t == "Camera":
            self.setCursor(Qt.ArrowCursor)
        elif t == "Select":
            self.setCursor(Qt.CrossCursor)
        elif t == "Modify":
            self.setCursor(Qt.OpenHandCursor)

    def mousePressEvent(self,e):
        self.setFocus()
        self.mousePos = e.pos()
        if self.mouseMode=='Modify':
            #initiate undo
            self.mol.initUndo()
            self.modMode=''
            #determine which atoms to modify
            sel = self.mol.getSelection()
            if sel:
                self.modData=[set(i[0] for i in sel)]
            else:
                self.modData = [range(self.mol.nat)]
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
                def selToCoord(sel):
                    return self.mol.getAtom(sel[0])[1]
                pick = self.pickAtom(e)
                #picked atom
                if pick:
                    self.modData+=[self.mol.getAtom(pick[0])[1]]
                #com of selection
                elif sel:
                    coords=[selToCoord(i) for i in sel]
                    self.modData+=[(np.max(coords,axis=0)+np.min(coords,axis=0))/2]
                #com of whole cell
                else:
                    self.modData+= [self.mol.getCenter(True)]
        elif e.buttons()&2:
            if self.mouseMode =='Camera':
                self.rMatrix.setToIdentity()
            elif self.mouseMode == 'Select':
                self.mol.delSelection()
                self.parent.updateMol()
            self.update()

    def mouseMoveEvent(self,e):
        delta=e.pos()-self.mousePos
        if (e.buttons() & 1):
            if self.mouseMode == 'Camera':
                #rotate camera
                tmp = QMatrix4x4()
                tmp.rotate(delta.x(),0,1,0)
                tmp.rotate(delta.y(),1,0,0)
                self.rMatrix = tmp*self.rMatrix
                self.update()
            elif self.mouseMode=='Select':
                #draw rectangle if over threshold
                if delta.manhattanLength()>5:
                    self.rectPos=e.pos()
                    self.update()
                #retain original mousePos
                return
            elif self.mouseMode=='Modify':
                #rotate selected atoms around center
                atoms=self.modData[0]
                angle= abs(delta.x()) + abs(delta.y())
                axis = delta.y()*self.modData[1][0]-delta.x()*self.modData[1][1]
                shift=self.modData[2]
                self.mol.rotate(atoms,angle,axis,shift)
                self.parent.updateMol()
                self.modMode='rotate'
        elif e.buttons()&2 and self.mouseMode=='Modify':
            #shift selection in z-direction (cam space)
            atoms=self.modData[0]
            vec = self.modData[1][2]*(delta.x()+delta.y())
            self.mol.shift(atoms,vec*0.1)
            self.parent.updateMol()
            self.modMode='shift'
        elif (e.buttons() & 4):
            if self.mouseMode=='Camera':
                #shift camera
                self.xsh += delta.x()/10.
                self.ysh -= delta.y()/10.
                self.update()
            elif self.mouseMode=='Modify':
                #shift selection in camera-plane
                atoms=self.modData[0]
                vec = self.modData[1][0]*delta.x()+self.modData[1][1]*delta.y()
                self.mol.shift(atoms,vec*0.1)
                self.parent.updateMol()
                self.modMode='shift'
        self.mousePos = e.pos()

    def mouseReleaseEvent(self,e):
        if self.mouseMode=='Modify' and self.modMode:
            self.mol.saveUndo(self.modMode)
            self.modMode=''
            self.parent.updateMol()
        elif self.mouseMode=='Select' and e.button()&1:
            if self.rectPos:
                for i in self.pickAtom(self.mousePos,self.rectPos):
                    self.mol.addSelection(i)
                self.rectPos = None
                self.parent.updateMol()
            else:
                pick = self.pickAtom(e)
                if pick:
                    self.mol.addSelection(pick)
                    self.parent.updateMol()

    def pickAtom(self,c1,c2=None):
        #render with selectionmode
        self.paintStuff(True)
        #Wait for everything to render,configure memory alignment
        glFlush()
        glFinish()
        glPixelStorei(GL_UNPACK_ALIGNMENT,1)
        #determine size of picked area
        if c2:
            x = min(c1.x(),c2.x())
            y = self.height()-1-max(c1.y(),c2.y())
            w = max(1,abs(c1.x()-c2.x()))
            h = max(1,abs(c1.y()-c2.y()))
            color = (GLubyte*(4*w*h))(0)
        else:
            x = c1.x()
            y = self.height()-1-c1.y()
            w = 1
            h = 1
            color = (GLubyte*4)(0)
        #Read pixel(s) from GPU
        glReadPixels(x,y,w,h,GL_RGBA,GL_UNSIGNED_BYTE,color)
        #determine unique colors
        cset=set()
        for i in range(0,len(color),4):
            cset.add(tuple(color[i:i+4]))
        #ignore background
        cset.discard((255,255,255,0))
        #identify atoms
        atoms=[]
        if cset:
            mult=self.mult[0]*self.mult[1]*self.mult[2]
            for i in cset:
                idx = i[0] + 256*i[1] + 65536*i[2]
                if idx<len(self.atomsVBO):
                    realid = idx//mult
                    off = idx%mult
                    zoff = off%self.mult[2]
                    yoff = (off//self.mult[2])%self.mult[1]
                    xoff = ((off//self.mult[2])//self.mult[1])%self.mult[0]
                    atoms.append((realid,(xoff,yoff,zoff)))
        if c2:
            return atoms
        elif atoms:
            return atoms[0]
        else:
            return None

    def wheelEvent(self,e):
        delta = e.delta()
        #zoom with vertical wheel
        if e.orientation() & 2:
                if delta < 0:
                        self.distance *= 1.1
                elif delta > 0:
                        self.distance *= 0.9
                self.update()
        e.accept()

    ##########################################
    # MODIFY RENDER-STATE
    ##########################################

    def alignView(self):
        if self.sender().text()=='+x':
            self.rMatrix = QMatrix4x4([0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1])
        elif self.sender().text()=='-x':
            self.rMatrix = QMatrix4x4([0,-1,0,0,0,0,1,0,-1,0,0,0,0,0,0,1])
        elif self.sender().text()=='+y':
            self.rMatrix = QMatrix4x4([-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1])
        elif self.sender().text()=='-y':
            self.rMatrix = QMatrix4x4([1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,1])
        elif self.sender().text()=='+z':
            self.rMatrix.setToIdentity()
        elif self.sender().text()=='-z':
            self.rMatrix = QMatrix4x4([-1,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,1])
        self.update()

    ##########################################
    # UPDATE RENDER-DATA
    ##########################################

    def setMol(self,mol,mult):
        #prepare atoms and bonds for drawing
        if not mol:return

        #save for interaction
        self.mol=mol
        self.mult=mult
        #local variables for convenience
        atoms = mol.getAtoms(fmt='bohr')
        pse = mol.pse
        vec = mol.getVec()*mol.getCellDim('bohr')
        center = mol.getCenter(config["Rotate around COM"])
        bonds = mol.getBonds(config['Bond cutoff factor'])
        sel = mol.getSelection()
        isoVal = mol.getIsoVal()
        volPlane = mol.getVolPlane()
        milPlane = mol.getMilPlane()

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
        r=np.float32(config['Bond radius'])
        n=np.float32(0)
        one=np.float32(1)
        for j in [0,1,2,3,4,5,6,7]:
            #check if pbc part is necessary
            if mult&j!=j: continue
            for i in bonds[j]:
                #get positions of atoms
                a = atoms[i[0]][1]+i[2][0]
                b = atoms[i[1]][1]+i[2][1]
                #save colors
                c1 = list(map(np.float32,pse[atoms[i[0]][0]]['col']))
                c2 = list(map(np.float32,pse[atoms[i[1]][0]]['col']))
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
                    c=np.float32(np.cos(theta))
                    s=np.float32(-np.sin(theta))
                    ic=np.float32(1.-c)
                for idx,k in enumerate(off):
                    if j>0 and edge[idx]&j!=0:
                        continue
                    self.bondPosVBO.append([l*(ic*ax[0]*ax[0]+c),l*(ic*ax[0]*ax[1]-s*ax[2]),l*(ic*ax[0]*ax[2]+s*ax[1]),n,
                                r*(ic*ax[0]*ax[1]+s*ax[2]),r*(ic*ax[1]*ax[1]+c),r*(ic*ax[1]*ax[2]-s*ax[0]),n,
                                r*(ic*ax[0]*ax[2]-s*ax[1]),r*(ic*ax[1]*ax[2]+s*ax[0]),r*(ic*ax[2]*ax[2]+c),n,
                                pos[0]+k[0],pos[1]+k[1],pos[2]+k[2],one,
                                c1[0],c1[1],c1[2],c1[3],
                                c2[0],c2[1],c2[2],c2[3]])
        if self.instanced:
            self.bondPosVBO=VBO(np.array(self.bondPosVBO,'f'))

        #save atoms in VBOs
        if config['Atom radius VdW']:
            rad = 'vdwr'
        else:
            rad = 'covr'
        self.atomsVBO=[(at[1]+j).tolist()+[pse[at[0]][rad],pse[at[0]]['col'][0],pse[at[0]]['col'][1],pse[at[0]]['col'][2],pse[at[0]]['col'][3]] for at in atoms for j in off]
        if self.instanced:
            self.atomsVBO=VBO(np.array(self.atomsVBO,'f'))
        #check for selected atoms inside mult-range
        if sel:
            self.selVBO=[]
            for i in sel:
                if all(i[1]<np.array(self.mult)):
                    at=atoms[i[0]]
                    pos = (at[1]+np.dot(i[1],vec)+off[0]).tolist()
                    self.selVBO.append(pos+[pse[at[0]][rad]*1.5,0.4,0.4,0.5,0.5])
            if self.selVBO and self.instanced:
                self.selVBO=VBO(np.array(self.selVBO,'f'))
        elif hasattr(self,'selVBO'):
            del self.selVBO
        #check for isoValue -> make isoSurface
        if isoVal[0]:
            self.surfVBO = VBO(makeIsoSurf(self.mol.getVol(),self.mol.getVolGradient(),isoVal[1],isoVal[2]))
        elif hasattr(self,'surfVBO'):
            del self.surfVBO
        #check for Volume-heatmap:
        if volPlane[0]:
            v = self.mol.getVol()
            vmin = v.min()
            vdiff = v.max()-vmin
            if volPlane[1]==0:
                pdat=v[volPlane[2],:,:]
                pos = volPlane[2]/v.shape[0]
                p=np.array([[pos,0,0],[pos,1,0],[pos,0,1],[pos,1,0],[pos,0,1],[pos,1,1]],'f')
            elif volPlane[1]==1:
                pdat=v[:,volPlane[2],:]
                pos = volPlane[2]/v.shape[1]
                p=np.array([[0,pos,0],[1,pos,0],[0,pos,1],[1,pos,0],[0,pos,1],[1,pos,1]],'f')
            elif volPlane[1]==2:
                pdat=v[:,:,volPlane[2]]
                pos = volPlane[2]/v.shape[2]
                p=np.array([[0,0,pos],[1,0,pos],[0,1,pos],[1,0,pos],[0,1,pos],[1,1,pos]],'f')
            self.volPlaneTex=np.array([(x-vmin)/vdiff for x in pdat],'f')
            vec=self.mol.getVec()*self.mol.getCellDim()
            UV = [[0,0],[0,1],[1,0],[0,1],[1,0],[1,1]]
            self.volPlaneVBO=VBO(np.array([np.dot(p[i],vec).tolist()+UV[i%6] for i in range(len(p))],'f'))
        elif hasattr(self,'volPlaneVBO'):
            del self.volPlaneTex
            del self.volPlaneVBO
        #check for crystal-plane:
        if milPlane[0]:
            self.milPlaneTex=np.array([[1.]],'f')
            p=[]
            #catch undefined case
            pval = milPlane[1]
            if pval.count(0) == 3:
                if hasattr(self,'milPlaneVBO'):
                    del self.milPlaneVBO
                    del self.milPlaneTex
                self.update()
                return
            elif pval[0] == 0:
                if pval[1] == 0:
                    for l in range(abs(pval[2])):
                        p+=[[0,0,(l+1.)/pval[2]],[1,0,(l+1.)/pval[2]],[0,1,(l+1.)/pval[2]],[1,0,(l+1.)/pval[2]],[0,1,(l+1.)/pval[2]],[1,1,(l+1.)/pval[2]]]
                elif pval[2] == 0:
                    for k in range(abs(pval[1])):
                        p+=[[0,(k+1.)/pval[1],0],[1,(k+1.)/pval[1],0],[0,(k+1.)/pval[1],1],[1,(k+1.)/pval[1],0],[0,(k+1.)/pval[1],1],[1,(k+1.)/pval[1],1]]
                else:
                    for k in range(abs(pval[1])):
                        for l in range(abs(pval[2])):
                            p+=[[0,float(k)/pval[1],(l+1.)/pval[2]],[0,(k+1.)/pval[1],float(l)/pval[2]],[1,float(k)/pval[1],(l+1.)/pval[2]],[0,(k+1.)/pval[1],float(l)/pval[2]],[1,float(k)/pval[1],(l+1.)/pval[2]],[1,(k+1.)/pval[1],float(l)/pval[2]]]
            else:
                if pval[1] == 0:
                    if pval[2]==0:
                        for h in range(abs(pval[0])):
                            p+=[[(h+1.)/pval[0],0,0],[(h+1.)/pval[0],1,0],[(h+1.)/pval[0],0,1],[(h+1.)/pval[0],1,0],[(h+1.)/pval[0],0,1],[(h+1.)/pval[0],1,1]]
                    else:
                        for h in range(abs(pval[0])):
                            for l in range(abs(pval[2])):
                                p+=[[(h+1.)/pval[0],0,float(l)/pval[2]],[(h+1.)/pval[0],1,float(l)/pval[2]],[float(h)/pval[0],0,(l+1.)/pval[2]],[(h+1.)/pval[0],1,float(l)/pval[2]],[float(h)/pval[0],0,(l+1.)/pval[2]],[float(h)/pval[0],1,(l+1.)/pval[2]]]
                elif pval[2] == 0:
                    for h in range(abs(pval[0])):
                        for k in range(abs(pval[1])):
                            p+=[[(h+1.)/pval[0],float(k)/pval[1],0],[(h+1.)/pval[0],float(k)/pval[1],1],[float(h)/pval[0],(k+1.)/pval[1],0],[(h+1.)/pval[0],float(k)/pval[1],1],[float(h)/pval[0],(k+1.)/pval[1],0],[float(h)/pval[0],(k+1.)/pval[1],1]]
                else:
                    for h in range(abs(pval[0])):
                        for k in range(abs(pval[1])):
                            for l in range(abs(pval[2])):
                                p+=[[(h+1.)/pval[0],float(k)/pval[1],float(l)/pval[2]],[float(h)/pval[0],(k+1.)/pval[1],float(l)/pval[2]],[float(h)/pval[0],float(k)/pval[1],(l+1.)/pval[2]],[(h+1.)/pval[0],(k+1.)/pval[1],float(l)/pval[2]],[(h+1.)/pval[0],float(k)/pval[1],(l+1.)/pval[2]],[float(h)/pval[0],(k+1.)/pval[1],(l+1.)/pval[2]]]
            p=np.array(p,'f')
            #take care of negative hkl-values
            for i in range(3):
                if pval[i]<0:
                    p[:,i]+=1

            #generate planeVBO
            vec=self.mol.getVec()*self.mol.getCellDim()
            UV = [[0,0],[0,1],[1,0],[0,1],[1,0],[1,1]]
            self.milPlaneVBO=VBO(np.array([np.dot(p[i],vec).tolist()+UV[i%6] for i in range(len(p))],'f'))
        elif hasattr(self,'milPlaneVBO'):
            del self.milPlaneVBO
            del self.milPlaneTex
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

        self.update()

    ###################################################
    # RENDERING
    ###################################################

    def initializeGL(self):
        #set line width for cell
        glLineWidth(2)
        glPointSize(2)

        #find supported version and modify shaders accordingly
        glVersion=float(glGetString(GL_VERSION)[0:3])
        self.instanced = glVersion>=3.3
        if self.instanced:
            self.glslv='#version 330\n'
            #prepare VAO
            self.VAO = glGenVertexArrays(1)
        else:
            self.glslv='#version 130\n'

        def makeShader(vf,ff):
            s = QGLShaderProgram()
            with open(dirname(__file__)+'/opengl/'+vf) as f:
                v=self.glslv+f.read()
                s.addShaderFromSourceCode(QGLShader.Vertex,v)
            with open(dirname(__file__)+'/opengl/'+ff) as f:
                f=self.glslv+f.read()
                s.addShaderFromSourceCode(QGLShader.Fragment,f)
            return s

        #add shaders:
        self.sphereShader = makeShader('vertexSpheres.vert','fragmentSpheres.frag')
        self.bondShader   = makeShader('vertexBonds.vert','fragmentBonds.frag')
        self.lineShader   = makeShader('vertexLines.vert','fragmentLines.frag')
        self.selectShader = makeShader('vertexSelect.vert','fragmentSelect.frag')
        self.planeShader  = makeShader('vertexPlane.vert','fragmentPlane.frag')
        self.surfShader   = makeShader('vertexSurf.vert','fragmentSpheres.frag')

        # load sphere
        sf=open(dirname(__file__)+'/opengl/sphere_model','r')
        self.sphereVBO = VBO(np.array(sf.readline().split(),'f'))
        sf.close()
        # load torus
        tf=open(dirname(__file__)+'/opengl/bond_model','r')
        self.torusVBO=VBO(np.array(tf.readline().split(),'f'))
        tf.close()

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

    def paintEvent(self,e):
        if self.rectPos:
            p = QPainter(self)
            if config['Antialiasing']:
                p.setRenderHint(QPainter.Antialiasing,True)
            self.paintStuff()
            #push transparent area in front
            glDisable(GL_DEPTH_TEST)
            glDisable(GL_CULL_FACE)
            #draw rect
            pen = QPen()
            pen.setWidth(2)
            p.setPen(pen)
            b = QBrush()
            b.setColor(QColor(180,180,180,40))
            b.setStyle(Qt.SolidPattern)
            p.setBrush(b)
            p.drawRect(QRectF(self.mousePos,self.rectPos))
            p.end()
        else:
            self.paintStuff()
            self.updateGL()

    def paintStuff(self,select=False):
        if self.instanced: glBindVertexArray(self.VAO)
        #clear depth and color buffer:
        self.qglClearColor(QColor(255,255,255,0))
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #set global GL settings
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_CULL_FACE)
        glEnable(GL_BLEND)
        if config['Antialiasing']:
            glEnable(GL_MULTISAMPLE)
        else:
            glDisable(GL_MULTISAMPLE)
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)

        if not hasattr(self,'atomsVBO'): return

        #construct viewMatrix
        self.vMatrix = QMatrix4x4()
        self.vMatrix.lookAt(QVector3D(0,0,self.distance),QVector3D(0,0,0),QVector3D(0,1,0))
        self.vMatrix.translate(self.xsh,self.ysh,0)

        #check for projection:
        if config["Perspective projection"]:
            self.proj = self.pMatrix
        else:
            self.proj = self.oMatrix
            #scale based on distance for zoom effect
            self.vMatrix.scale(10./self.distance)
        #rendering:
        if select:
            glDisable(GL_MULTISAMPLE)
            self.drawAtomsSelect()
        else:
            self.drawAtoms()
            if config["Show bonds"]:
                self.drawBonds()
            if config["Show cell"]:
                self.drawCell()
            if hasattr(self,'surfVBO'):
                self.drawSurf()
            if hasattr(self,'volPlaneVBO'):
                self.drawPlane(self.volPlaneVBO,self.volPlaneTex)
            if hasattr(self,'milPlaneVBO'):
                self.drawPlane(self.milPlaneVBO,self.milPlaneTex)
            if hasattr(self,'selVBO'):
                self.drawSelection()
        if self.instanced: glBindVertexArray(0)

    def drawAtoms(self):
        self.sphereShader.bind()

        self.sphereVBO.bind()
        glEnableVertexAttribArray(0)
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,None)
        self.sphereVBO.unbind()

        self.sphereShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
        self.sphereShader.setUniformValue('rMatrix',self.rMatrix)
        self.sphereShader.setUniformValue('atom_fac',config["Atom radius factor"])

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
            glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO)//3,len(self.atomsVBO))
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
                glDrawArrays(GL_TRIANGLES,0,len(self.sphereVBO)//3)

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

        self.surfShader.setUniformValue('volOff',*self.mol.getVolOffset().tolist())
        self.surfShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)
        self.surfShader.setUniformValue('cellVec',QMatrix3x3((self.mol.getVec()*self.mol.getCellDim()).flatten()))
        self.surfShader.setUniformValue('rMatrix',self.rMatrix)

        glDisable(GL_CULL_FACE)

        if self.instanced:
            self.offVBO.bind()
            glEnableVertexAttribArray(3)
            glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,0,None)
            glVertexAttribDivisor(3,1)
            self.offVBO.unbind()
            glDrawArraysInstanced(GL_TRIANGLES,0,len(self.surfVBO)//9,len(self.offVBO))
            glVertexAttribDivisor(3,0)
            glDisableVertexAttribArray(3)
        else:
            for i in self.offVBO:
                self.surfShader.setUniformValue('offset',*i)
                glDrawArrays(GL_TRIANGLES,0,len(self.surfVBO)//9)

        glEnable(GL_CULL_FACE)
        glDisableVertexAttribArray(0)
        glDisableVertexAttribArray(1)
        glDisableVertexAttribArray(2)
        self.surfShader.release()

    def drawPlane(self,plane,tex):
        self.planeShader.bind()

        plane.bind()
        glEnableVertexAttribArray(0)
        glEnableVertexAttribArray(1)
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,20,None)
        glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,20,plane+12)
        plane.unbind()

        self.planeShader.setUniformValue('vpMatrix',self.proj*self.vMatrix*self.rMatrix)

        #send texture
        ID = glGenTextures(1)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D,ID)
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
        glTexImage2D(GL_TEXTURE_2D,0,GL_RED,tex.shape[1],tex.shape[0],0,GL_RED,GL_FLOAT,tex)
        self.planeShader.setUniformValue('texSampler',0)

        glDisable(GL_CULL_FACE)

        if self.instanced:
            self.offVBO.bind()
            glEnableVertexAttribArray(2)
            glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,0,None)
            glVertexAttribDivisor(2,1)
            self.offVBO.unbind()
            glDrawArraysInstanced(GL_TRIANGLES,0,len(plane),len(self.offVBO))
            glVertexAttribDivisor(2,0)
            glDisableVertexAttribArray(2)
        else:
            for i in self.offVBO:
                self.planeShader.setUniformValue('offset',*i)
                glDrawArrays(GL_TRIANGLES,0,len(plane))

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
        self.selectShader.setUniformValue('atom_fac',config["Atom radius factor"])

        if self.instanced:
            self.atomsVBO.bind()
            glEnableVertexAttribArray(1)
            glEnableVertexAttribArray(2)
            glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,None)
            glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,self.atomsVBO+12)
            glVertexAttribDivisor(1,1)
            glVertexAttribDivisor(2,1)
            self.atomsVBO.unbind()
            glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO)//3,len(self.atomsVBO))
            glDisableVertexAttribArray(1)
            glDisableVertexAttribArray(2)
            glVertexAttribDivisor(1,0)
            glVertexAttribDivisor(2,0)
        else:
            for j,i in enumerate(self.atomsVBO):
                self.selectShader.setUniformValue('position_modelspace',*i[0:3])
                self.selectShader.setUniformValue('scale_modelspace',i[3])
                self.selectShader.setUniformValue('in_color',(j&0xFF)/255.,((j&0xFF00)>>8)/255.,((j&0xFF0000)>>16)/255.,1)
                glDrawArrays(GL_TRIANGLES,0,len(self.sphereVBO)//3)

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
        self.sphereShader.setUniformValue('atom_fac',config["Atom radius factor"])

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
            glDrawArraysInstanced(GL_TRIANGLES,0,len(self.sphereVBO)//3,len(self.selVBO))
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
                glDrawArrays(GL_TRIANGLES,0,len(self.sphereVBO)//3)

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
            glDrawArraysInstanced(GL_TRIANGLES,0,len(self.torusVBO)//3,len(self.bondPosVBO))
            for i in range(1,7):
                glDisableVertexAttribArray(i)
                glVertexAttribDivisor(i,0)
        else:
            for i in self.bondPosVBO:
                self.bondShader.setUniformValue('mMatrix',QMatrix4x4(i[:16]).transposed())
                self.bondShader.setUniformValue('s1Color',*i[16:20])
                self.bondShader.setUniformValue('s2Color',*i[20:])
                glDrawArrays(GL_TRIANGLES,0,len(self.torusVBO)//3)

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
