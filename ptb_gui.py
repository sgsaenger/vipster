#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

from copy import deepcopy
from os.path import splitext
from os import getcwd
#import numpy as np

from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer,Qt

from viewport import ViewPort
from coordedit import MolArea
from paramedit import PWTab
from multiedit import ToolArea

class MainWindow(QMainWindow):

        def __init__(self,controller):
                super(MainWindow,self).__init__()
                self.controller = controller
                self.initApp()

        def initApp(self):

                #Create Menu
                self.initMenu()

                #Create main widget
                mv = MainView(self,self.controller)
                self.setCentralWidget(mv)

                # Set Title and run:
                self.setWindowTitle('PWToolBox')
                self.show()

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
                exitAction = QAction('&Exit',self)
                exitAction.setShortcut('Ctrl+Q')
                exitAction.triggered.connect(qApp.quit)

                fMenu = self.menuBar().addMenu('&File')
                fMenu.addAction(newAction)
                fMenu.addAction(loadAction)
                fMenu.addAction(saveAction)
                fMenu.addSeparator()
                fMenu.addAction(exitAction)

	def newHandler(self):
		self.controller.newMol()
		self.centralWidget().loadView()

        def loadHandler(self):
                fname = QFileDialog.getOpenFileName(self,'Open File',getcwd())
                if not fname: return
                ftype = QInputDialog.getItem(self,'Choose file type','File type:',self.controller.indict.keys(),0,False)
                ftype = str(ftype[0])
                self.controller.readFile(ftype,fname)
                self.centralWidget().loadView()

        def saveHandler(self):
                fname = QFileDialog.getSaveFileName(self,'Save File',getcwd())
                if not fname: return
                ftype = QInputDialog.getItem(self,'Choose File type','File type:',self.controller.outdict.keys(),0,False)
                ftype = str(ftype[0])
                mol = self.centralWidget().getMolecule()
                if ftype=='PWScf Input':
                        param = self.centralWidget().getParam()
                else:
                        param = False
                coordfmt = self.centralWidget().coord.fmt.currentText()
                self.controller.writeFile(ftype,mol,fname,param,coordfmt)

class MainView(QWidget):

        def __init__(self,parent,controller):
                super(MainView,self).__init__()
                self.parent = parent
                self.controller = controller
                # initialize GUI and accompanying actions
                self.initUI()
                self.initActions()
                # set persistent render parameters
                self.mult =[1,1,1]
                self.view = True
                self.cell = True
                self.mouseMode = True
                self.AA = True
                self.bondShow = True
                self.selection=[]

        def initUI(self):

        #Left column:
                #Molecule list:
                self.mlist = QListWidget()
                self.mlist.currentRowChanged.connect(self.selectMolecule)
                #splitter-container
                mlist=QWidget()
                mlayout=QVBoxLayout()
                mlabel=QLabel('Loaded Molecules:')
                mlayout.addWidget(mlabel)
                mlayout.addWidget(self.mlist)
                mlist.setLayout(mlayout)

                #PWParameter list:
                self.pwlist = QListWidget()
                self.pwlist.currentRowChanged.connect(self.selectPWParam)
                #splitter-container
                pwlist=QWidget()
                pwlayout=QVBoxLayout()
                pwlabel=QLabel('PW Parameter sets:')
                pwlayout.addWidget(pwlabel)
                pwlayout.addWidget(self.pwlist)
                pwlist.setLayout(pwlayout)


                #Edit stuff
                self.edit = ToolArea(self)

                #encapsulate in splitter:
                lcol = QSplitter()
                lcol.addWidget(mlist)
                lcol.addWidget(pwlist)
                lcol.addWidget(self.edit)
                lcol.setOrientation(0)
                lcol.setChildrenCollapsible(False)
                lcol.setFrameStyle(38)
                lcol.setMaximumWidth(300)


        #Central Column:
                #OpenGL Viewport:
                self.visual = ViewPort(self)

                #Control visual:
                #Cell multiplication
                self.xspin = QSpinBox()
                self.yspin = QSpinBox()
                self.zspin = QSpinBox()
                for i in [self.xspin,self.yspin,self.zspin]:
                        i.setMinimum(1)
                #Toggle edit column:
                editBut = QPushButton()
                editBut.setText('Edit')
                editBut.setCheckable(True)
                #create Timers
                self.animTimer = QTimer()
                self.animTimer.setInterval(50)
                self.animTimer.timeout.connect(self.incStep)
                #Screenshot test
                screenbut = QPushButton()
                screenbut.setText('Save Screenshot')
                screenbut.clicked.connect(self.makeScreen)
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
                self.currentStep = QLineEdit()
                self.currentStep.setText('0')
                self.currentStep.setValidator(QIntValidator(0,0))
                self.maxStep = QLabel('0')
                #connect updateMolStep as everything is ready
                self.currentStep.textChanged.connect(self.updateMolStep)
                self.xspin.valueChanged.connect(self.updateMult)
                self.yspin.valueChanged.connect(self.updateMult)
                self.zspin.valueChanged.connect(self.updateMult)
                #Control Layout:
                mult = QHBoxLayout()
                mult.addWidget(QLabel('Cell multiply:'))
                mult.addStretch()
                mult.addWidget(QLabel('x:'))
                mult.addWidget(self.xspin)
                mult.addStretch()
                mult.addWidget(QLabel('y:'))
                mult.addWidget(self.yspin)
                mult.addStretch()
                mult.addWidget(QLabel('z:'))
                mult.addWidget(self.zspin)
                mult.addWidget(screenbut)
                steps = QHBoxLayout()
                steps.addWidget(firstBut)
                steps.addWidget(decBut)
                steps.addWidget(self.currentStep)
                steps.addWidget(QLabel('/'))
                steps.addWidget(self.maxStep)
                steps.addWidget(playBut)
                steps.addWidget(incBut)
                steps.addWidget(lastBut)
                steps.addWidget(editBut)

                #Layout:
                viewlay = QVBoxLayout()
                viewlay.addWidget(self.visual)
                viewlay.addLayout(mult)
                viewlay.addLayout(steps)

                #Frame it:
                mcol = QFrame()
                mcol.setLayout(viewlay)
                mcol.setFrameStyle(38)

        #Right column:
                #Molecule edit area:
                self.coord = MolArea(self)

                #PWParameter edit area:
                self.pw = PWTab()

                #nest edit areas in tabwidget
                rcol = QTabWidget()
                rcol.addTab(self.coord,'Molecule coordinates')
                rcol.addTab(self.pw,'PW Parameters')
                rcol.setFixedWidth(467)
                #connect to toggle button
                rcol.hide()
                editBut.clicked.connect(rcol.setVisible)

        #Lay out columns:
                hbox = QHBoxLayout()
                hbox.addWidget(lcol)
                hbox.addWidget(mcol)
                hbox.addWidget(rcol)
                self.setLayout(hbox)

        def initActions(self):
                #define actions
                changeProj = QAction('Change &Projection',self)
                changeProj.setShortcut('p')
                changeProj.triggered.connect(self.projHandler)
                toggleCell = QAction('Toggle &Cell',self)
                toggleCell.setShortcut('c')
                toggleCell.triggered.connect(self.cellHandler)
                mouseRotate = QAction('&Rotate',self)
                mouseRotate.setShortcut('r')
                mouseRotate.triggered.connect(self.mmHandler)
                mouseSelect = QAction('&Select',self)
                mouseSelect.setShortcut('s')
                mouseSelect.triggered.connect(self.mmHandler)
                #antiAlias = QAction('Toggle &Antialiasing',self)
                #antiAlias.setShortcut('a')
                #antiAlias.triggered.connect(self.aaHandler)
                toggleBonds = QAction('Toggle &Bonds',self)
                toggleBonds.setShortcut('b')
                toggleBonds.triggered.connect(self.bondHandler)
                #create menu
                vMenu = self.parent.menuBar().addMenu('&View')
                vMenu.addAction(changeProj)
                vMenu.addAction(toggleCell)
                vMenu.addAction(toggleBonds)
                #vMenu.addAction(antiAlias)
                vMenu.addSeparator()
                vMenu.addAction(mouseRotate)
                vMenu.addAction(mouseSelect)

        ########################################################
        #insert loaded molecules
        ########################################################
        def loadView(self):
                count = self.mlist.count()
                for i in range(count,self.controller.get_nmol()):
                        self.mlist.addItem("Mol "+str(i+1))
                self.mlist.setCurrentRow(self.mlist.count()-1)
                count = self.pwlist.count()
                for i in range(count,self.controller.get_npw()):
                        self.pwlist.addItem("Param "+str(i+1))
                self.pwlist.setCurrentRow(self.pwlist.count()-1)

        ########################################################
        #get data from controller
        ########################################################
        def selectMolecule(self,sel):
                steps = self.controller.get_lmol(sel)
                self.currentStep.setValidator(QIntValidator(1,steps))
                self.currentStep.setText(str(steps))
                self.maxStep.setText(str(steps))
                self.updateMolStep()

        def selectPWParam(self,sel):
                self.pw.setPW(self.controller.get_pw(sel))

        def updateMolStep(self):
                #get current Molecule from MolList
                sel = self.mlist.currentRow()
                step = int(self.currentStep.text())-1
                mol = self.controller.get_mol(sel,step)
                #Send Molecule to Visualisation and Editor
                self.coord.setMol(mol)
                self.edit.setMol(mol)
                self.visual.setMol(mol,self.mult)

        def updateMult(self):
                self.mult=[self.xspin.value(),self.yspin.value(),self.zspin.value()]
                self.updateMolStep()

        ########################################################
        #to controller
        ########################################################
        def getMolecule(self):
                return self.controller.get_mol(self.mlist.currentRow(),int(self.currentStep.text())-1)

        def getParam(self):
                self.pw.saveParam()
                return self.controller.get_pw(self.pwlist.currentRow())

        ########################################################
        # view modifier
        ########################################################

        def projHandler(self):
                self.view = not self.view
                self.visual.updateGL()

        def cellHandler(self):
                self.cell = not self.cell
                self.visual.updateGL()

        ########################################################
        # Toggle Bonds
        ########################################################

        def bondHandler(self):
                self.bondShow = not self.bondShow
                self.visual.updateGL()

        ########################################################
        # Mouse mode
        ########################################################

        def mmHandler(self):
                source = self.sender()
                if source == None: self.mouseMode = not self.mouseMode
                elif source.text() == '&Rotate': self.mouseMode = True
                elif source.text() == '&Select': self.mouseMode = False

        ########################################################
        # Toggle AntiAliasing
        ########################################################

        def aaHandler(self):
                if self.AA:
                        self.AA = False
                        glDisable(GL_MULTISAMPLE)
                elif not self.AA:
                        self.AA = True
                        glEnable(GL_MULTISAMPLE)
                self.visual.updateGL()

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
                if self.currentStep.text()==self.maxStep.text():
                        self.animTimer.stop()
                else:
                        self.currentStep.setText(str(int(self.currentStep.text())+1))

        def decStep(self):
                if self.currentStep.text()=='1':return
                self.currentStep.setText(str(int(self.currentStep.text())-1))

        def firstStep(self):
                self.currentStep.setText('1')

        def lastStep(self):
                self.currentStep.setText(self.maxStep.text())

        def toggleAnim(self):
                if self.animTimer.isActive():
                        self.animTimer.stop()
                else:
                        self.animTimer.start()

