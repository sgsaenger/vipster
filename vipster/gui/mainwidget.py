# -*- coding: utf-8 -*-

from os.path import splitext
from os import getcwd

from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer,Qt

from .viewport import ViewPort
from .moltab import MolTab
from .pwtab import PWTab
from .lammpstab import LammpsTab
from .conftab import ConfTab
from .tools import tools
from .collapsiblewidget import collapsibleWidget

from .. import *
from ..ftypeplugins import _gui_indict,_gui_outdict

GuiMolecules = []
GuiParameters = []

def launchVipster(m=GuiMolecules,p=GuiParameters):
    app = QApplication([])
    gui = MainWidget(m,p)
    app.aboutToQuit.connect(gui.deleteLater)
    app.exec_()

__all__=['GuiMolecules','GuiParameters','launchVipster']

class MainWidget(QWidget):

        def __init__(self,m,p):
            super(MainWidget,self).__init__()
            self.parent = QMainWindow()
            self.molecules=m
            self.parameters=p
            self.initMenu()
            self.mult =[1,1,1]

        #Right column:
            rcol = QTabWidget()
            #Molecule edit area:
            self.moltab = MolTab(self)
            rcol.addTab(self.moltab,"Molecule coordinates")
            #Parameter area
            self.pwtab = PWTab()
            self.lammpstab = LammpsTab()
            self.noparam=QLabel("No Parameter set selected")
            paramWidget=QWidget()
            paramLayout=QVBoxLayout()
            paramLayout.addWidget(self.noparam)
            paramLayout.addWidget(self.pwtab,stretch=1)
            paramLayout.addWidget(self.lammpstab,stretch=0)
            paramLayout.addStretch(0)
            self.pwtab.hide()
            self.lammpstab.hide()
            paramWidget.setLayout(paramLayout)
            rcol.addTab(paramWidget,"Parameters")
            #Settings:
            self.conftab = ConfTab(self)
            rcol.addTab(self.conftab,"Settings")
            #connect to toggle button
            rcol.setFixedWidth(467)
            rcol.hide()

        #Central Column:
            #OpenGL Viewport:
            self.visual = ViewPort(self)
            #mouse mode selection
            self.mouse = QButtonGroup()
            mouse = QHBoxLayout()
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
                mouse.addWidget(i)
                self.mouse.addButton(i)
                i.setCheckable(True)
            self.mouse.buttonClicked.connect(self.visual.setMouseMode)
            self.mouse.setExclusive(True)
            camBut.setChecked(True)
            #Cell multiplication
            self.xspin = QSpinBox()
            self.yspin = QSpinBox()
            self.zspin = QSpinBox()
            for i in [self.xspin,self.yspin,self.zspin]:
                    i.setMinimum(1)
            #Toggle right column:
            editBut = QPushButton()
            editBut.setText("Edit")
            editBut.setCheckable(True)
            editBut.clicked.connect(rcol.setVisible)
            #create Timers
            self.animTimer = QTimer()
            self.animTimer.setInterval(50)
            self.animTimer.timeout.connect(self.incStep)
            #camera-buttons
            xcam = QPushButton("+x")
            xcam.clicked.connect(self.visual.alignView)
            xcam.setFixedWidth(30)
            ycam = QPushButton("+y")
            ycam.clicked.connect(self.visual.alignView)
            ycam.setFixedWidth(30)
            zcam = QPushButton("+z")
            zcam.clicked.connect(self.visual.alignView)
            zcam.setFixedWidth(30)
            mxcam = QPushButton("-x")
            mxcam.clicked.connect(self.visual.alignView)
            mxcam.setFixedWidth(30)
            mycam = QPushButton("-y")
            mycam.clicked.connect(self.visual.alignView)
            mycam.setFixedWidth(30)
            mzcam = QPushButton("-z")
            mzcam.clicked.connect(self.visual.alignView)
            mzcam.setFixedWidth(30)
            #Choose step, animate
            incBut = QPushButton()
            incBut.clicked.connect(self.incStep)
            incBut.setAutoRepeat(True)
            incBut.setIcon(self.style().standardIcon(65))
            decBut = QPushButton()
            decBut.clicked.connect(self.decStep)
            decBut.setAutoRepeat(True)
            decBut.setIcon(self.style().standardIcon(66))
            playBut = QPushButton()
            playBut.setIcon(self.style().standardIcon(60))
            playBut.clicked.connect(self.toggleAnim)
            firstBut = QPushButton()
            firstBut.setIcon(self.style().standardIcon(64))
            firstBut.clicked.connect(self.firstStep)
            lastBut = QPushButton()
            lastBut.setIcon(self.style().standardIcon(63))
            lastBut.clicked.connect(self.lastStep)
            self.Step = QSlider()
            self.Step.setOrientation(1)
            self.Step.setMinimum(1)
            self.Step.setMaximum(1)
            self.Step.setTickPosition(self.Step.TicksBelow)
            self.Step.setSingleStep(1)
            self.Step.valueChanged.connect(self.updateMolStep)
            self.curStep = QLabel("1")
            self.maxStep = QLabel("1")
            self.xspin.valueChanged.connect(self.updateMult)
            self.yspin.valueChanged.connect(self.updateMult)
            self.zspin.valueChanged.connect(self.updateMult)
            #Control Layout:
            mult = QHBoxLayout()
            mult.addWidget(QLabel("Align camera:"))
            mult.addStretch(1)
            mult.addWidget(xcam)
            mult.addStretch(1)
            mult.addWidget(mxcam)
            mult.addStretch(1)
            mult.addWidget(ycam)
            mult.addStretch(1)
            mult.addWidget(mycam)
            mult.addStretch(1)
            mult.addWidget(zcam)
            mult.addStretch(1)
            mult.addWidget(mzcam)
            mult.addStretch(100)
            mult.addWidget(QLabel("Cell multiply:"))
            mult.addStretch(1)
            mult.addWidget(QLabel("x:"))
            mult.addWidget(self.xspin)
            mult.addStretch(1)
            mult.addWidget(QLabel("y:"))
            mult.addWidget(self.yspin)
            mult.addStretch(1)
            mult.addWidget(QLabel("z:"))
            mult.addWidget(self.zspin)
            steps = QHBoxLayout()
            steps.addWidget(firstBut)
            steps.addWidget(decBut)
            steps.addWidget(self.Step)
            steps.addWidget(self.curStep)
            steps.addWidget(QLabel("/"))
            steps.addWidget(self.maxStep)
            steps.addWidget(playBut)
            steps.addWidget(incBut)
            steps.addWidget(lastBut)
            steps.addWidget(editBut)
            #Layout:
            viewlay = QVBoxLayout()
            viewlay.addLayout(mouse)
            viewlay.addWidget(self.visual)
            viewlay.addLayout(mult)
            viewlay.addLayout(steps)
            #Frame it:
            mcol = QFrame()
            mcol.setLayout(viewlay)
            mcol.setFrameStyle(38)

        #Left column:
            #Molecule list:
            self.mlist = QListWidget()
            self.mlist.currentRowChanged.connect(self.selectMolecule)
            self.mlist.setMinimumHeight(150)
            for i in self.molecules:
                self.mlist.addItem(i.name)
            #PWParameter list:
            self.paramlist = QListWidget()
            self.paramlist.currentRowChanged.connect(self.selectParam)
            self.paramlist.setMinimumHeight(150)
            for i in self.parameters:
                self.paramlist.addItem(i["name"])
            #Layout and frame
            lcol = QScrollArea()
            lcol.setFrameStyle(38)
            lcol.setFixedWidth(340)
            llay = QVBoxLayout()
            llay.addWidget(QLabel("Loaded Molecules:"))
            llay.addWidget(self.mlist)
            llay.addWidget(collapsibleWidget("Parameter sets:",self.paramlist,show=False))
            self.edit=[]
            for i in tools.items():
                self.edit.append(i[1](self))
                llay.addWidget(collapsibleWidget(i[0],self.edit[-1],show=False))
            lcolarea= QWidget()
            lcolarea.setLayout(llay)
            lcolarea.setFixedWidth(320)
            lcol.setWidget(lcolarea)
            lcol.setWidgetResizable(True)
            lcol.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        #Lay out columns:
            hbox = QHBoxLayout()
            hbox.addWidget(lcol)
            hbox.addWidget(mcol)
            hbox.addWidget(rcol)
            self.setLayout(hbox)

        #Show window
            self.parent.setCentralWidget(self)
            self.parent.setWindowTitle("PWToolBox")
            self.parent.show()
            self.mlist.setCurrentRow(self.mlist.count()-1)
            self.paramlist.setCurrentRow(self.paramlist.count()-1)

        def initMenu(self):
            fMenu = self.parent.menuBar().addMenu("&File")
            newAction = QAction("&New Molecule",self)
            newAction.setShortcut("Ctrl+N")
            newAction.triggered.connect(self.newMolHandler)
            fMenu.addAction(newAction)
            pMenu = fMenu.addMenu("New &Parameter set")
            newPW = QAction("PWScf Parameter set",self)
            newPW.triggered.connect(self.newPWHandler)
            pMenu.addAction(newPW)
            newLammps = QAction("LAMMPS Parameter set",self)
            newLammps.triggered.connect(self.newLammpsHandler)
            pMenu.addAction(newLammps)
            loadAction = QAction("&Load Molecule(s)",self)
            loadAction.setShortcut("Ctrl+O")
            loadAction.triggered.connect(self.loadHandler)
            fMenu.addAction(loadAction)
            saveAction = QAction("&Save Molecule",self)
            saveAction.setShortcut("Ctrl+S")
            saveAction.triggered.connect(self.saveHandler)
            fMenu.addAction(saveAction)
            scrotAction = QAction("Save Screensho&t",self)
            scrotAction.setShortcut("Ctrl+P")
            scrotAction.triggered.connect(self.makeScreen)
            fMenu.addAction(scrotAction)
            fMenu.addSeparator()
            exitAction = QAction("&Exit",self)
            exitAction.setShortcut("Ctrl+Q")
            exitAction.triggered.connect(qApp.quit)
            fMenu.addAction(exitAction)

        def newMolHandler(self):
            self.molecules.append(Molecule())
            self.mlist.addItem("New Mol")
            self.mlist.setCurrentRow(self.mlist.count()-1)

        def newPWHandler(self):
            param = newParam("pwi")
            self.parameters.append(param)
            self.paramlist.addItem(param["name"])
            self.paramlist.setCurrentRow(self.paramlist.count()-1)

        def newLammpsHandler(self):
            param = newParam("lmp")
            self.parameters.append(param)
            self.paramlist.addItem(param["name"])
            self.paramlist.setCurrentRow(self.paramlist.count()-1)

        def loadHandler(self):
            fname = QFileDialog.getOpenFileName(self,"Open File",getcwd())
            if not fname: return
            ftype = QInputDialog.getItem(self,"Choose file type","File type:",list(_gui_indict.keys()),0,False)
            if not ftype[1]: return
            ftype = str(ftype[0])
            m,p = readFile(fname,ftype,mode="gui")
            self.molecules.append(m)
            self.mlist.addItem("Mol")
            self.mlist.setCurrentRow(self.mlist.count()-1)
            if p:
                self.parameters.append(p)
                self.paramlist.addItem(p["name"])
                self.paramlist.setCurrentRow(self.paramlist.count()-1)

        def saveHandler(self):
            fname = QFileDialog.getSaveFileName(self,"Save File",getcwd())
            if not fname: return
            ftype = QInputDialog.getItem(self,"Choose File type","File type:",list(_gui_outdict.keys()),0,False)
            if not ftype[1]: return
            ftype = str(ftype[0])
            try:
                mol = self.curMol
                if ftype=="PWScf Input":
                    try:
                        param = self.parameters[self.paramlist.currentRow()]
                        if param["type"]!="pw":raise IndexError
                    except:
                        raise IndexError("No PW Parameter set")
                elif ftype=="Lammps Data file":
                    try:
                        param = self.parameters[self.paramlist.currentRow()]
                        if param["type"]!="lammps":raise IndexError
                    except:
                        raise IndexError("No LAMMPS Parameter set")
                else:
                    param = False
                coordfmt = self.moltab.fmt.currentText()
                writeFile(mol,ftype,fname,param,coordfmt,mode="gui")
            except StandardError as e:
                QMessageBox(QMessageBox.Critical,"Error",type(e).__name__+": "+e.message,QMessageBox.Ok,self).exec_()

        ########################################################
        #update view upon selections
        ########################################################
        def selectMolecule(self,sel):
            self.curMol = self.molecules[sel]
            steps=len(self.curMol)
            self.maxStep.setText(str(steps))
            self.Step.setMaximum(steps)
            self.Step.setValue(steps)
            self.updateMolStep()

        def selectParam(self,sel):
            param = self.parameters[sel]
            if param["type"]=="pw":
                self.noparam.hide()
                self.lammpstab.hide()
                self.pwtab.show()
                self.pwtab.setParam(param)
            elif param["type"]=="lammps":
                self.noparam.hide()
                self.pwtab.hide()
                self.lammpstab.show()
                self.lammpstab.setParam(param)

        def updateMolStep(self):
            #change step of trajectory when needed
            step = self.Step.value()-1
            self.curStep.setText(str(step+1))
            self.curMol.changeStep(step)
            #Send Molecule to Visualisation and Editor
            self.moltab.setMol(self.curMol)
            self.visual.setMol(self.curMol,self.mult)
            self.conftab.setMol(self.curMol)
            for i in self.edit:
                i.setMol(self.curMol)

        def updateMult(self):
            self.mult=[self.xspin.value(),self.yspin.value(),self.zspin.value()]
            self.updateMolStep()

        ########################################################
        #screenshot test
        ########################################################
        def makeScreen(self):
                img = self.visual.grabFrameBuffer(True)
                fn = QFileDialog.getSaveFileName(self,"Save Screenshot",getcwd(),"Portable Network Graphics (*.png)")
                if not fn: return
                if splitext(str(fn))[1] == "": fn+=".png"
                img.save(fn,None,0)

        ########################################################
        #steps and animation:
        ########################################################
        def incStep(self):
                if self.Step.value()==self.Step.maximum():
                        self.animTimer.stop()
                else:
                        self.Step.setValue(self.Step.value()+1)

        def decStep(self):
                if self.Step.value()==1:return
                self.Step.setValue(self.Step.value()-1)

        def firstStep(self):
                self.Step.setValue(1)

        def lastStep(self):
                self.Step.setValue(self.Step.maximum())

        def toggleAnim(self):
                if self.animTimer.isActive():
                        self.animTimer.stop()
                else:
                        self.animTimer.start()

