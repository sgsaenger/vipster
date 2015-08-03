#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

from os.path import splitext
from os import getcwd

from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer

from viewport import ViewPort
from moltab import MolTab
from pwtab import PWTab
from tools import tools
from collapsiblewidget import collapsibleWidget

class MainView(QWidget):

        def __init__(self,controller):
            super(MainView,self).__init__()
            self.parent = QMainWindow()
            self.controller = controller
            self.initMenu()
            self.mult =[1,1,1]

        #Right column:
            #Molecule edit area:
            self.coord = MolTab(self)
            #PWParameter edit area:
            self.pw = PWTab()
            #nest edit areas in tabwidget
            rcol = QTabWidget()
            rcol.addTab(self.coord,'Molecule coordinates')
            rcol.addTab(self.pw,'PW Parameters')
            rcol.setFixedWidth(467)
            #connect to toggle button
            rcol.hide()

        #Central Column:
            #OpenGL Viewport:
            self.visual = ViewPort(self)
            #mouse mode selection
            self.mouse = QButtonGroup()
            mouse = QHBoxLayout()
            camBut = QPushButton('Camera')
            camBut.setShortcut('r')
            camBut.setToolTip('Move camera\n\nLMB: Rotate camera\nMMB: Drag camera\nRMB: Align camera to z-axis')
            selBut = QPushButton('Select')
            selBut.setShortcut('s')
            selBut.setToolTip('Select atoms\n\nLMB: Select atoms\nRMB: Clear selection')
            modBut = QPushButton('Modify')
            modBut.setShortcut('m')
            modBut.setToolTip('Modify geometry\n\nLMB: Rotate atoms (around Center of Mass or selected Atom)\nMMB: Move atoms in xy-plane (camera)\nRMB: Move atoms along z-axis (camera)')
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
            editBut.setText('Edit')
            editBut.setCheckable(True)
            editBut.clicked.connect(rcol.setVisible)
            #create Timers
            self.animTimer = QTimer()
            self.animTimer.setInterval(50)
            self.animTimer.timeout.connect(self.incStep)
            #camera-buttons
            xcam = QPushButton('x')
            xcam.clicked.connect(self.visual.alignView)
            xcam.setFixedWidth(30)
            ycam = QPushButton('y')
            ycam.clicked.connect(self.visual.alignView)
            ycam.setFixedWidth(30)
            zcam = QPushButton('z')
            zcam.clicked.connect(self.visual.alignView)
            zcam.setFixedWidth(30)
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
            self.curStep = QLabel('1')
            self.maxStep = QLabel('1')
            self.xspin.valueChanged.connect(self.updateMult)
            self.yspin.valueChanged.connect(self.updateMult)
            self.zspin.valueChanged.connect(self.updateMult)
            #Control Layout:
            mult = QHBoxLayout()
            mult.addWidget(QLabel('Align camera:'))
            mult.addStretch(1)
            mult.addWidget(xcam)
            mult.addStretch(1)
            mult.addWidget(ycam)
            mult.addStretch(1)
            mult.addWidget(zcam)
            mult.addStretch(100)
            mult.addWidget(QLabel('Cell multiply:'))
            mult.addStretch(1)
            mult.addWidget(QLabel('x:'))
            mult.addWidget(self.xspin)
            mult.addStretch(1)
            mult.addWidget(QLabel('y:'))
            mult.addWidget(self.yspin)
            mult.addStretch(1)
            mult.addWidget(QLabel('z:'))
            mult.addWidget(self.zspin)
            steps = QHBoxLayout()
            steps.addWidget(firstBut)
            steps.addWidget(decBut)
            steps.addWidget(self.Step)
            steps.addWidget(self.curStep)
            steps.addWidget(QLabel('/'))
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
            #PWParameter list:
            self.pwlist = QListWidget()
            self.pwlist.currentRowChanged.connect(self.selectPWParam)
            #Layout and frame
            lcol = QFrame()
            lcol.setFrameStyle(38)
            lcol.setFixedWidth(340)
            llay = QVBoxLayout()
            llay.addWidget(QLabel('Loaded Molecules:'))
            llay.addWidget(self.mlist)
            llay.addWidget(collapsibleWidget('PW Parameter sets:',self.pwlist,show=False))
            self.edit=[]
            for i in tools.items():
                self.edit.append(i[1](self))
                llay.addWidget(collapsibleWidget(i[0],self.edit[-1],show=False))
            lcol.setLayout(llay)

        #Lay out columns:
            hbox = QHBoxLayout()
            hbox.addWidget(lcol)
            hbox.addWidget(mcol)
            hbox.addWidget(rcol)
            self.setLayout(hbox)

        #Show window
            self.parent.setCentralWidget(self)
            self.parent.setWindowTitle('PWToolBox')
            self.parent.show()

        def initMenu(self):
            newAction = QAction('&New Molecule',self)
            newAction.setShortcut('Ctrl+N')
            newAction.triggered.connect(self.newHandler)
            loadAction = QAction('&Load Molecule(s)',self)
            loadAction.setShortcut('Ctrl+O')
            loadAction.triggered.connect(self.loadHandler)
            saveAction = QAction('&Save Molecule',self)
            saveAction.setShortcut('Ctrl+S')
            saveAction.triggered.connect(self.saveHandler)
            scrotAction = QAction('Save Screensho&t',self)
            scrotAction.setShortcut('Ctrl+P')
            scrotAction.triggered.connect(self.makeScreen)
            exitAction = QAction('&Exit',self)
            exitAction.setShortcut('Ctrl+Q')
            exitAction.triggered.connect(qApp.quit)

            fMenu = self.parent.menuBar().addMenu('&File')
            fMenu.addAction(newAction)
            fMenu.addAction(loadAction)
            fMenu.addAction(saveAction)
            fMenu.addAction(scrotAction)
            fMenu.addSeparator()
            fMenu.addAction(exitAction)

        def newHandler(self):
            self.controller.newMol()
            self.loadView()

        def loadHandler(self):
            fname = QFileDialog.getOpenFileName(self,'Open File',getcwd())
            if not fname: return
            ftype = QInputDialog.getItem(self,'Choose file type','File type:',self.controller.gui_indict.keys(),0,False)
            ftype = str(ftype[0])
            self.controller.readFile(ftype,fname)
            self.loadView()

        def saveHandler(self):
            fname = QFileDialog.getSaveFileName(self,'Save File',getcwd())
            if not fname: return
            ftype = QInputDialog.getItem(self,'Choose File type','File type:',self.controller.gui_outdict.keys(),0,False)
            ftype = str(ftype[0])
            try:
                mol = self.getMolecule()
                if ftype=='PWScf Input':
                    try:
                        param = self.getParam()
                    except:
                        raise IndexError('No PW Parameter set')
                else:
                        param = False
                coordfmt = self.coord.fmt.currentText()
                self.controller.writeFile(ftype,mol,fname,param,coordfmt)
            except StandardError as e:
                QMessageBox(QMessageBox.Critical,'Error',e.message,QMessageBox.Ok,self).exec_()

        ########################################################
        #insert loaded molecules
        ########################################################
        def loadView(self):
                count = self.mlist.count()
                for i in range(count,self.controller.getNMol()):
                        self.mlist.addItem("Mol "+str(i+1))
                self.mlist.setCurrentRow(self.mlist.count()-1)
                count = self.pwlist.count()
                for i in range(count,self.controller.getNPw()):
                        self.pwlist.addItem("Param "+str(i+1))
                self.pwlist.setCurrentRow(self.pwlist.count()-1)

        ########################################################
        #get data from controller
        ########################################################
        def selectMolecule(self,sel):
            self.mol = self.controller.getMol(sel)
            steps=len(self.mol)
            self.maxStep.setText(str(steps))
            self.Step.setMaximum(steps)
            self.Step.setValue(steps)
            self.updateMolStep()

        def selectPWParam(self,sel):
            self.pw.setPW(self.controller.getPw(sel))

        def updateMolStep(self):
            #change step of trajectory when needed
            step = self.Step.value()-1
            self.curStep.setText(str(step+1))
            self.mol.changeMol(step)
            #Send Molecule to Visualisation and Editor
            self.coord.setMol(self.mol)
            self.visual.setMol(self.mol,self.mult)
            for i in self.edit:
                i.setMol(self.mol)

        def updateMult(self):
            self.mult=[self.xspin.value(),self.yspin.value(),self.zspin.value()]
            self.updateMolStep()

        ########################################################
        #to controller
        ########################################################
        def getMolecule(self):
                return self.controller.getMol(self.mlist.currentRow())

        def getParam(self):
                return self.controller.getPw(self.pwlist.currentRow())


        ########################################################
        #screenshot test
        ########################################################
        def makeScreen(self):
                img = self.visual.grabFrameBuffer(True)
                fn = QFileDialog.getSaveFileName(self,'Save Screenshot',getcwd(),'Portable Network Graphics (*.png)')
                if not fn: return
                if splitext(str(fn))[1] == '': fn+='.png'
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

