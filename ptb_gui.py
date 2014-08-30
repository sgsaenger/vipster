#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

from copy import deepcopy
from os.path import dirname,splitext
from os import getcwd
from math import floor,sqrt
import numpy as np

from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer,Qt
from PyQt4.QtOpenGL import *
from OpenGL.GL import *

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
                self.initActions()
                menu = self.menuBar()
                fMenu = menu.addMenu('&File')
                fMenu.addAction(self.loadAction)
                fMenu.addAction(self.saveAction)
                fMenu.addSeparator()
                fMenu.addAction(self.exitAction)

        def initActions(self):
                self.loadAction = QAction('&Load',self)
                self.loadAction.setShortcut('Ctrl+O')
                self.loadAction.triggered.connect(self.loadHandler)
                self.saveAction = QAction('&Save',self)
                self.saveAction.setShortcut('Ctrl+S')
                self.saveAction.triggered.connect(self.saveHandler)
                self.exitAction = QAction('&Exit',self)
                self.exitAction.setShortcut('Ctrl+Q')
                self.exitAction.triggered.connect(qApp.quit)

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
                param = self.centralWidget().getParam()
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
                # set persisten render parameters
                self.mult =[1,1,1]
                self.view = True
                self.cell = True
                self.mouseMode = True
                self.selection=[]

        def initUI(self):

        #Left column:
                #Molecule list:
                self.mlist = QListWidget()
                self.mlist.currentRowChanged.connect(self.updateMolList)
                #splitter-container
                mlist=QWidget()
                mlayout=QVBoxLayout()
                mlabel=QLabel('Loaded Molecules:')
                mlayout.addWidget(mlabel)
                mlayout.addWidget(self.mlist)
                mlist.setLayout(mlayout)

                #PWParameter list:
                self.pwlist = QListWidget()
                self.pwlist.currentRowChanged.connect(self.setParam)
                #splitter-container
                pwlist=QWidget()
                pwlayout=QVBoxLayout()
                pwlabel=QLabel('PW Parameter sets:')
                pwlayout.addWidget(pwlabel)
                pwlayout.addWidget(self.pwlist)
                pwlist.setLayout(pwlayout)


                #Edit stuff
                self.edit = EditArea(self)

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
                #create menu
                vMenu = self.parent.menuBar().addMenu('&View')
                vMenu.addAction(changeProj)
                vMenu.addAction(toggleCell)
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
        def updateMolList(self,sel):
                steps = self.controller.get_lmol(sel)
                self.currentStep.setValidator(QIntValidator(1,steps))
                self.currentStep.setText(str(steps))
                self.maxStep.setText(str(steps))
                self.updateMolStep()

        def setParam(self,sel):
                self.pw.setPW(self.controller.get_pw(sel))

        def updateMolStep(self):
                #get current Molecule from MolList
                sel = self.mlist.currentRow()
                step = int(self.currentStep.text())-1
                mol = self.controller.get_mol(sel,step)
                #Send Molecule to Visualisation and Editor
                self.coord.setMol(mol)
                self.edit.setMol(mol)
                self.setMol(mol)

        def updateMult(self):
                self.setMult([self.xspin.value(),self.yspin.value(),self.zspin.value()])

        #################################################
        # CALLED UPON SELECTING MOLECULE
        #################################################

        def setMol(self,mol):
                self.mol = mol
                self.makeCell()
                self.prepObjects()

        def makeCell(self):
                #get vectors:
                vec = self.mol.get_vec()
                cdm = self.mol.get_celldm()
                null = QVector3D(0,0,0)
                a = QVector3D(vec[0,0],vec[0,1],vec[0,2])*cdm
                b = QVector3D(vec[1,0],vec[1,1],vec[1,2])*cdm
                c = QVector3D(vec[2,0],vec[2,1],vec[2,2])*cdm
                #define vertices:
                self.cell_modelspace=[null,a,null,b,null,c,
                                      a,a+b,a,a+c,
                                      b,b+a,b,b+c,
                                      c,c+a,c,c+b,
                                      a+b,a+b+c,
                                      a+c,a+b+c,
                                      b+c,a+b+c]

        def setMult(self,mult):
                self.mult = mult
                self.prepObjects()

        def prepObjects(self):
                if not hasattr(self,'mol'):return
                #local variables for convenience
                mult = self.mult
                atoms = [self.mol.get_atom(i) for i in range(self.mol.get_nat())]
                vec = self.mol.get_vec()*self.mol.get_celldm()
                center = self.mol.get_center()

                #get bonds and offsets
                if mult == [1,1,1]:
                        #self.off = [-self.mol.get_center()]
                        self.off = [[0,0,0]]
                        bonds = [self.mol.get_bonds(),[],[],[],[],[],[],[]]
                else:
                        self.off = []
                        tmult = [1,1,1]
                        #save the multiplicators:
                        for i in [0,1,2]:
                                if self.mult[i]%2 == 0:
                                        tmult[i]=[x+0.5-self.mult[i]/2 for x in range(self.mult[i])]
                                else:
                                        tmult[i]=[x-floor(self.mult[i]/2) for x in range(self.mult[i])]
                        #generate offsets:
                        for i in tmult[0]:
                                for j in tmult[1]:
                                        for k in tmult[2]:
                                                self.off.append([i,j,k])
                        #self.off = self.getOffsets()
                        bonds = self.mol.get_pbc_bonds()

                #prepare bonds with position,angle,axis and names (for coloring)
                tempbonds=[[],[],[],[],[],[],[],[]]
                for j in range(8):
                        for i in bonds[j]:
                                #get positions of atoms
                                a = i[0]
                                b = i[1]
                                #save colors
                                c1 = self.visual.pse[i[2]][2]
                                c2 = self.visual.pse[i[3]][2]
                                #position of bond
                                pos= (a+b)/2
                                #rotate bond from x-axis d to bond-axis c
                                c = (a-b)/np.linalg.norm(a-b)
                                d = np.array([1,0,0])
                                theta = np.degrees(np.arccos(np.dot(c,d)))
                                axis = -np.cross(c,d)
                                tempbonds[j].append([pos,theta,axis,c1,c2])

                #save positions of atoms and bonds in world space
                self.atoms=[]
                self.bonds = []
                self.cells = []
                edge = np.max(self.off,axis=0)
                for i1 in self.off:
                        off = (i1[0]*vec[0]+i1[1]*vec[1]+i1[2]*vec[2])-center
                        self.cells.append(off)
                        for j in atoms:
                                #save coord,color,size
                                self.atoms.append((j[1]+off,self.visual.pse[j[0]][2],self.visual.pse[j[0]][1]))
                        for j in tempbonds[0]:
                                self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                                pass
                        #vec0
                        if mult[0]!=1 and i1[0]!=edge[0]:
                                for j in tempbonds[1]:
                                        self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                                #vec0+vec1
                                if mult[1]!=1 and i1[1]!=edge[1]:
                                        for j in tempbonds[4]:
                                                self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                                        #vec0+vec1+vec2
                                        if mult[2]!=1 and i1[2]!=edge[2]:
                                                for j in tempbonds[7]:
                                                        self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                                #vec0+vec2
                                if mult[2]!=1 and i1[2]!=edge[2]:
                                        for j in tempbonds[5]:
                                                self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                        #vec1
                        if mult[1]!=1 and i1[1]!=edge[1]:
                                for j in tempbonds[2]:
                                        self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                                #vec1+vec2
                                if mult[2]!=1 and i1[2]!=edge[2]:
                                        for j in tempbonds[6]:
                                                self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                        #vec2
                        if mult[2]!=1 and i1[2]!=edge[2]:
                                for j in tempbonds[3]:
                                        self.bonds.append([j[0]+off,j[1],j[2],j[3],j[4]])
                self.visual.updateGL()

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
        # Mouse mode
        ########################################################

        def mmHandler(self):
                source = self.sender()
                if source == None: self.mouseMode = not self.mouseMode
                elif source.text() == '&Rotate': self.mouseMode = True
                elif source.text() == '&Select': self.mouseMode = False

        ########################################################
        #screenshot test
        ########################################################
        def makeScreen(self):
                img = self.visual.grabFrameBuffer(True)
                fn = QFileDialog.getSaveFileName(self,'Save Screenshot',getcwd(),'Portable Network Graphics (*.png)')
                if not fn: return
                if splitext(str(fn))[1] == '': fn+='.png'
                img.save(fn,None,100)

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

class MolArea(QWidget):

        def __init__(self,parent):
                super(MolArea,self).__init__()
                self.initTab()
                self.parent = parent

        def initTab(self):
                # coord fmt dropdown selector
                self.fmt = QComboBox()
                self.fmt.setToolTip('Select format of coordinates')
                for i in ['angstrom','bohr','crystal','alat']:
                        self.fmt.addItem(i)
                self.fmt.setCurrentIndex(2)
                self.fmt.currentIndexChanged.connect(self.fillTab)

                # layout1
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('Coordinates:'))
                hbox.addWidget(self.fmt)

                # show coordinates in table
                self.table = QTableWidget()
                self.table.setColumnCount(4)
                self.table.setHorizontalHeaderLabels(['Type','x','y','z'])
                self.table.itemChanged.connect(self.cellHandler)
                #show actions in right click menu
                self.table.setContextMenuPolicy(Qt.ActionsContextMenu)

                # show celldm
                self.cellDm = QLineEdit()
                self.cellDm.editingFinished.connect(self.updateCellDm)

                # show cell vectors in table
                self.vtable = QTableWidget()
                self.vtable.setColumnCount(3)
                self.vtable.setRowCount(3)
                self.vtable.setFixedHeight(120)
                self.vtable.setHorizontalHeaderLabels(['x','y','z'])
                self.vtable.itemChanged.connect(self.vecHandler)

                #cell header with label and celldm
                hbox2=QHBoxLayout()
                hbox2.addWidget(QLabel('Cell vectors:'))
                hbox2.addStretch(1)
                hbox2.addWidget(QLabel('Cell dimension:'))
                hbox2.addWidget(self.cellDm)

                #New Atom button:
                newB = QPushButton('New Atom')
                newB.clicked.connect(self.newAtom)

                #Copy Atom button:
                copyB = QPushButton('Copy Atom(s)')
                copyB.clicked.connect(self.copyAt)

                #paste button:
                pasteB = QPushButton('Paste Atom(s)')
                pasteB.clicked.connect(self.pasteAt)

                #delete button:
                delB = QPushButton('Delete Atom(s)')
                delB.clicked.connect(self.delAt)

                #sort buttons
                btns = QHBoxLayout()
                btns.addWidget(copyB)
                btns.addWidget(pasteB)
                btns.addWidget(newB)
                btns.addWidget(delB)

                # actions
                self.newA = QAction('New Atom',self)
                self.newA.setShortcut('Ctrl+N')
                self.newA.triggered.connect(self.newAtom)
                self.table.addAction(self.newA)
                self.copyA = QAction('Copy Atom(s)',self)
                self.copyA.setShortcut('Ctrl+C')
                self.copyA.triggered.connect(self.copyAt)
                self.table.addAction(self.copyA)
                self.pasteA = QAction('Paste Atom(s)',self)
                self.pasteA.setShortcut('Ctrl+V')
                self.pasteA.triggered.connect(self.copyAt)
                self.table.addAction(self.pasteA)
                self.delA = QAction('Delete Atom(s)',self)
                self.delA.setShortcut('Del')
                self.delA.triggered.connect(self.delAt)
                self.table.addAction(self.delA)

                # Action modifiers
                #TODO: cellHandler,newAt,delAt,pasteAt ?
                self.appAll = QCheckBox('Apply to all Molecules')
                self.scale = QCheckBox('Scale coordinates with cell')
                hbox3 = QHBoxLayout()
                hbox3.addWidget(self.appAll)
                hbox3.addWidget(self.scale)

                # set Layout for Tab
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(self.table)
                vbox.addLayout(btns)
                vbox.addLayout(hbox2)
                vbox.addWidget(self.vtable)
                vbox.addLayout(hbox3)

                self.setLayout(vbox)
                self.resize(self.sizeHint())

        #load selected molecule
        def setMol(self,mol):
                #connect molecule
                self.mol = mol
                #fill if visible
                if self.isVisible(): self.fillTab()

        def showEvent(self,e):
                if hasattr(self,'mol'): self.fillTab()

        ##################################################################
        # EDIT FUNCTIONS
        ##################################################################
        def newAtom(self):
                self.mol.create_atom()
                self.mol.set_bonds()
                self.mol.set_pbc_bonds()

                #update Main Widget
                self.parent.updateMolStep()

        def copyAt(self):
                self.sel = []
                for i in self.table.selectedRanges():
                        for j in range(i.topRow(),i.bottomRow()+1):
                                self.sel.append(self.mol.get_atom(j))

        def pasteAt(self):
                pos = self.table.currentRow()+1
                self.sel.reverse()
                for at in self.sel:
                        self.mol.insert_atom(pos,at)
                self.mol.set_bonds()
                self.mol.set_pbc_bonds()

                #update Main Widget
                self.parent.updateMolStep()

        def delAt(self):
                delrange = set()
                for i in self.table.selectedRanges():
                        for j in range(i.topRow(),i.bottomRow()+1):
                                delrange.add(j)
                for i in sorted(delrange,reverse=True):
                        self.mol.del_atom(i)
                self.mol.set_bonds()
                self.mol.set_pbc_bonds()
                #update Main Widget
                self.parent.updateMolStep()

        #####################################################
        # UPDATE HANDLER
        ####################################################

        def updateCellDm(self):
                if self.updatedisable: return
                if float(self.cellDm.text()) == self.mol.get_celldm():return
                if self.appAll.isChecked():
                        if self.scale.isChecked():
                                par=self.parent
                                for i in range(int(par.maxStep.text())):
                                        mol=par.controller.get_mol(par.mlist.currentRow(),i)
                                        molc = deepcopy(mol)
                                        mol.set_celldm(float(self.cellDm.text()))
                                        for j in range(mol.get_nat()):
                                                mol.set_atom(j,*molc.get_atom(j,'crystal'))
                                        mol.set_bonds()
                                        mol.set_pbc_bonds()
                                        del(molc)
                        else:
                                par = self.parent
                                for i in range(int(par.maxStep.text())):
                                        par.controller.get_mol(par.mlist.currentRow(),i).set_celldm(float(self.cellDm.text()))
                else:
                        if self.scale.isChecked():
                                mol = deepcopy(self.mol)
                                self.mol.set_celldm(float(self.cellDm.text()))
                                for i in range(mol.get_nat()):
                                        self.mol.set_atom(i,*mol.get_atom(i,'crystal'))
                                self.mol.set_bonds()
                                self.mol.set_pbc_bonds()
                                del(mol)
                        else:
                            self.mol.set_celldm(float(self.cellDm.text()))

                #update Main Widget
                self.parent.updateMolStep()

        def cellHandler(self):
                if self.updatedisable: return
                atom = self.table.currentRow()
                name = str(self.table.item(atom,0).text())
                coord = [0,0,0]
                for j in [0,1,2]:
                        coord[j]=float(self.table.item(atom,j+1).text())
                self.mol.set_atom(atom,name,coord,self.fmt.currentText())
                self.mol.set_bonds()
                self.mol.set_pbc_bonds()

                #update Main Widget
                self.parent.updateMolStep()

        def vecHandler(self):
                if self.updatedisable: return
                vec=[[0,0,0],[0,0,0],[0,0,0]]
                for i in [0,1,2]:
                        for j in [0,1,2]:
                                vec[i][j]=float(self.vtable.item(i,j).text())
                if vec == self.mol.get_vec().tolist(): return
                if self.appAll.isChecked():
                        if self.scale.isChecked():
                                par = self.parent
                                for i in range(int(par.maxStep.text())):
                                        mol = par.controller.get_mol(par.mlist.currentRow(),i)
                                        molc = deepcopy(mol)
                                        mol.set_vec(vec)
                                        for j in range(mol.get_nat()):
                                                mol.set_atom(j,*molc.get_atom(j,'crystal'))
                                        mol.set_bonds()
                                        mol.set_pbc_bonds()
                                        del(molc)
                        else:
                                par = self.parent
                                for i in range(int(par.maxStep.text())):
                                        par.controller.get_mol(par.mlist.currentRow(),i).set_vec(vec)
                else:
                        if self.scale.isChecked():
                                mol = deepcopy(self.mol)
                                self.mol.set_vec(vec)
                                for i in range(mol.get_nat()):
                                        self.mol.set_atom(i,*mol.get_atom(i,'crystal'))
                                self.mol.set_bonds()
                                self.mol.set_pbc_bonds()
                                del(mol)
                        else:
                                self.mol.set_vec(vec)

                #update Main Widget
                self.parent.updateMolStep()

        ##############################################################
        # MAIN WIDGET UPDATE FUNCTION
        #############################################################

        def fillTab(self):
                #prevent handling of cell changes during fill
                self.updatedisable = True
                self.cellDm.setText(str(self.mol.get_celldm()))
                #fill atom table
                self.table.setRowCount(self.mol.get_nat())
                for i in range(self.mol.get_nat()):
                        at = self.mol.get_atom(i,self.fmt.currentText())
                        self.table.setItem(i,0,QTableWidgetItem(at[0]))
                        for j in [0,1,2]:
                                self.table.setItem(i,j+1,QTableWidgetItem(str(at[1][j])))
                #fill cell vec list
                vec = self.mol.get_vec()
                for i in [0,1,2]:
                        for j in [0,1,2]:
                                self.vtable.setItem(i,j,QTableWidgetItem(str(vec[i,j])))
                #reenable handling
                self.updatedisable = False

class PWTab(QSplitter):

        def __init__(self):
                super(PWTab,self).__init__()
                self.initTab()

        def initTab(self):
                self.tree = self.Tree()
                self.initKpoints()
                self.setOrientation(0)
                self.setChildrenCollapsible(False)
                self.setFixedWidth(445)
                self.addWidget(self.tree)
                self.addWidget(self.kp)
                self.setStretchFactor(0,1)

        #call fill functions when necessary:
        def setPW(self,pw):
                #save previous changes if applicable:
                if self.isVisible():
                        self.saveParam()
                #load parameters
                self.pw = pw
                #show if visible
                self.fillTree()
                self.fillKpoints()

        def showEvent(self,e):
                if hasattr(self,'pw'):
                        self.fillTree()
                        self.fillKpoints()

        def hideEvent(self,e):
                self.saveParam()

        def setCurrentIndex(self,i):
                if i>2: i=2
                self.kp.disp.setCurrentIndex(i)

        #init sections:
        def initKpoints(self):
                self.kp = QWidget()
                kp = self.kp
                #choose k point format
                kp.fmt = QComboBox()
                for i in ['gamma','automatic','tpiba','crystal','tpiba_b','crystal_b']:
                        kp.fmt.addItem(i)

                #Automatic: x,y,z, offset(x,y,z)
                kp.auto = QWidget()
                kp.auto.widg = []
                for i in [0,1,2]:
                        kp.auto.widg.append(QLineEdit())
                        kp.auto.widg[-1].setValidator(QIntValidator(1,40000))
                for i in [0,1,2]:
                        kp.auto.widg.append(QCheckBox())
                kp.auto.hbox1=QHBoxLayout()
                kp.auto.hbox1.addWidget(QLabel('x:'))
                kp.auto.hbox1.addWidget(kp.auto.widg[0])
                kp.auto.hbox1.addWidget(QLabel('y:'))
                kp.auto.hbox1.addWidget(kp.auto.widg[1])
                kp.auto.hbox1.addWidget(QLabel('z:'))
                kp.auto.hbox1.addWidget(kp.auto.widg[2])
                kp.auto.hbox2 = QHBoxLayout()
                kp.auto.hbox2.addWidget(QLabel('x offs.:'))
                kp.auto.hbox2.addWidget(kp.auto.widg[3])
                kp.auto.hbox2.addWidget(QLabel('y offs.:'))
                kp.auto.hbox2.addWidget(kp.auto.widg[4])
                kp.auto.hbox2.addWidget(QLabel('z offs.:'))
                kp.auto.hbox2.addWidget(kp.auto.widg[5])
                kp.auto.vbox = QVBoxLayout()
                kp.auto.vbox.addLayout(kp.auto.hbox1)
                kp.auto.vbox.addLayout(kp.auto.hbox2)
                kp.auto.setLayout(kp.auto.vbox)

                #discrete kpoints
                kp.disc = QTableWidget()
                kp.disc.setColumnCount(4)
                kp.disc.setHorizontalHeaderLabels(['x','y','z','weight'])
                kp.disc.setContextMenuPolicy(Qt.ActionsContextMenu)
                newKp = QAction('New k-point',kp.disc)
                newKp.setShortcut('Ctrl+K')
                newKp.triggered.connect(self.newKpoint)
                kp.disc.addAction(newKp)
                delKp = QAction('Delete k-point',kp.disc)
                delKp.triggered.connect(self.delKpoint)
                kp.disc.addAction(delKp)

                #stacked display of various formats
                kp.disp = QStackedWidget()
                kp.disp.addWidget(QLabel('Gamma point only'))
                kp.disp.addWidget(kp.auto)
                kp.disp.addWidget(kp.disc)
                kp.fmt.currentIndexChanged.connect(self.setCurrentIndex)

                #layout
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('K Points:'))
                hbox.addWidget(kp.fmt)
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(kp.disp)
                kp.setLayout(vbox)

        class Tree(QTreeWidget):
                def __init__(self):
                        super(PWTab.Tree,self).__init__()
                        self.setColumnCount(2)
                        self.setHeaderLabels(['Parameter','Value'])
                        self.setContextMenuPolicy(Qt.ActionsContextMenu)

                        # Actions:
                        newNl = QAction('New Namelist',self)
                        newNl.setShortcut('Ctrl+N')
                        newNl.triggered.connect(self.createNamelist)
                        self.addAction(newNl)
                        newPar = QAction('New Parameter',self)
                        newPar.setShortcut('Ctrl+P')
                        newPar.triggered.connect(self.createParameter)
                        self.addAction(newPar)
                        delItem = QAction('Delete Item',self)
                        delItem.setShortcut('Del')
                        delItem.triggered.connect(self.deleteItem)
                        self.addAction(delItem)

                def mouseDoubleClickEvent(self,e):
                        #on double left click, edit selected item
                        if (e.buttons() & 1):
                                self.editItem(self.currentItem(),self.currentColumn())

                def createNamelist(self):
                        new = QTreeWidgetItem(self)
                        new.setText(0,'&')
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                def createParameter(self):
                        #new = QTreeWidgetItem(self)
                        if self.currentItem().parent():
                                new = QTreeWidgetItem(self.currentItem().parent())
                        else:
                                new = QTreeWidgetItem(self.currentItem())
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                def deleteItem(self):
                        if self.currentItem().parent():
                                self.currentItem().parent().removeChild(self.currentItem())
                        else:
                                self.invisibleRootItem().removeChild(self.currentItem())

        #fill sections:
        def fillTree(self):
                root = self.tree.invisibleRootItem()
                #delete previous entries
                for i in range(root.childCount()):
                        root.removeChild(root.child(0))
                #mandatory namelists
                for i in ['&control','&system','&electrons']:
                        new = QTreeWidgetItem(self.tree)
                        new.setText(0,i)
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                #show optional namelists only if existing
                if '&ions' in self.pw:
                        new = QTreeWidgetItem(self.tree)
                        new.setText(0,'&ions')
                        new.setFlags(new.flags()|Qt.ItemIsEditable)
                if '&cell' in self.pw:
                        new = QTreeWidgetItem(self.tree)
                        new.setText(0,'&cell')
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                #show child entries
                for i in range(root.childCount()):
                        for j in self.pw[str(root.child(i).text(0))].items():
                                new = QTreeWidgetItem(root.child(i))
                                new.setText(0,j[0])
                                new.setText(1,j[1])
                                new.setFlags(new.flags()|Qt.ItemIsEditable)
                self.tree.expandAll()
                self.tree.resizeColumnToContents(0)

        def fillKpoints(self):
                self.kp.fmt.setCurrentIndex(['gamma','automatic','tpiba','crystal','tpiba_b','crystal_b'].index(self.pw['K_POINTS']['active']))
                if 'automatic' in self.pw['K_POINTS']:
                        for i in [0,1,2]:
                                self.kp.auto.widg[i].setText(self.pw['K_POINTS']['automatic'][i])
                        for i in [3,4,5]:
                                self.kp.auto.widg[i].setChecked(bool(int(self.pw['K_POINTS']['automatic'][i])))
                else:
                        for i in [0,1,2]:
                                self.kp.auto.widg[i].setText('')
                        for i in [3,4,5]:
                                self.kp.auto.widg[i].setChecked(False)
                if 'disc' in self.pw['K_POINTS']:
                        self.kp.disc.setRowCount(len(self.pw['K_POINTS']['disc']))
                        for i in range(len(self.pw['K_POINTS']['disc'])):
                                for j in [0,1,2,3]:
                                        self.kp.disc.setItem(i,j,QTableWidgetItem(self.pw['K_POINTS']['disc'][i][j]))
                else:
                        self.kp.disc.setRowCount(0)

        def newKpoint(self):
                self.kp.disc.setRowCount(self.kp.disc.rowCount()+1)

        def delKpoint(self):
                tab = self.kp.disc
                tab.removeRow(tab.currentRow())

        def saveParam(self):
                if not hasattr(self,'pw'):return
                #save monkhorst-pack-grid
                auto=[str(self.kp.auto.widg[i].text()) for i in [0,1,2]]
                auto+=[int(self.kp.auto.widg[i].isChecked()) for i in [3,4,5]]
                self.pw['K_POINTS']['automatic']=auto
                #save discrete k-points
                if self.kp.disc.rowCount() > 0:
                        disc = []
                        for i in range(self.kp.disc.rowCount()):
                                disc.append([self.kp.disc.item(i,j).text() for j in [0,1,2,3]])
                        self.pw['K_POINTS']['disc']=disc
                #save chosen k-point format
                self.pw['K_POINTS']['active']=str(self.kp.fmt.currentText())
                #save NameLists and parameters:
                for i in range(self.tree.invisibleRootItem().childCount()):
                        nl = self.tree.invisibleRootItem().child(i)
                        self.pw[str(nl.text(0))]={}
                        for i in range(nl.childCount()):
                                self.pw[str(nl.text(0))][str(nl.child(i).text(0))]=str(nl.child(i).text(1))

class EditArea(QWidget):
        def __init__(self,parent):
                super(EditArea,self).__init__()
                self.initStack()
                self.initMult()
                self.parent = parent

        def setMol(self,mol):
                self.mol = mol

        def initStack(self):
                self.stack = QStackedWidget()
                self.combo = QComboBox()
                self.combo.currentIndexChanged.connect(self.stack.setCurrentIndex)
                vbox = QVBoxLayout()
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('Edit:'))
                hbox.addWidget(self.combo)
                vbox.addLayout(hbox)
                vbox.addWidget(self.stack)
                self.setLayout(vbox)

        def initMult(self):
                self.mult = QWidget()
                self.combo.addItem('Mult. cell')
                self.stack.addWidget(self.mult)
                self.xmult = QLineEdit()
                self.xmult.setText('1')
                self.xmult.setValidator(QIntValidator(0,999))
                self.ymult = QLineEdit()
                self.ymult.setText('1')
                self.ymult.setValidator(QIntValidator(0,999))
                self.zmult = QLineEdit()
                self.zmult.setText('1')
                self.zmult.setValidator(QIntValidator(0,999))
                multBut = QPushButton('Apply')
                multBut.clicked.connect(self.multHandler)
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('x:'))
                hbox.addWidget(self.xmult)
                hbox.addWidget(QLabel('y:'))
                hbox.addWidget(self.ymult)
                hbox.addWidget(QLabel('z:'))
                hbox.addWidget(self.zmult)
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(multBut)
                self.mult.setLayout(vbox)

        def multHandler(self):
                if not hasattr(self,'mol'): return
                self.mol.mult(int(self.xmult.text()),int(self.ymult.text()),int(self.zmult.text()))
                self.mol.set_bonds()
                self.mol.set_pbc_bonds()
                self.parent.updateMolStep()

class ViewPort(QGLWidget):

        ##################################################
        # CALLED UPON INITIALISATION
        ##################################################
        def __init__(self,parent,aa=False):
                #init with antialiasing
                super(ViewPort,self).__init__(QGLFormat(QGL.SampleBuffers|QGL.AlphaChannel))
                self.parent = parent
                self.sphereShader = QGLShaderProgram()
                self.lineShader = QGLShaderProgram()
                self.bondShader = QGLShaderProgram()
                self.selectShader = QGLShaderProgram()
                self.xsh = 0
                self.ysh = 0
                self.rotMat = QMatrix4x4()
                self.distance = 25
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
                #render only visible vertices
                glEnable(GL_DEPTH_TEST)
                #backface culling: render only front of vertex
                glEnable(GL_CULL_FACE)
                #enable transparency for selection
                glEnable(GL_BLEND)
                glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)

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

                # prepare sphere
                self.makeSphere()
                #prepare bond
                self.makeBond()

        ##########################
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
                for i in range(3):
                        sphere = self.sphere_subdiv(sphere)
                self.atom_modelspace = sphere

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
                for i in range(2):
                        cylinder = self.cylinder_subdiv(cylinder)
                cylinder = self.cylinder_puzzle(cylinder)
                self.bond_modelspace = cylinder
                #set normals
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
                if not hasattr(self.parent,'mol'):return

                #initialize transformation matrices:
                self.mMatrix = QMatrix4x4()
                self.vMatrix = QMatrix4x4()

                #camera stuff:
                cRot = QMatrix4x4()
                cTrans = QMatrix4x4()

                #perspective
                self.vMatrix.lookAt(QVector3D(0,0,self.distance),QVector3D(0,0,0),QVector3D(0,1,0))

                #translate the view point:
                self.vMatrix.translate(self.xsh,self.ysh,0)

                #apply rotation
                self.vMatrix*=self.rotMat

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
                if select:
                        self.drawAtomsSelect()
                else:
                        self.drawAtoms()
                        self.drawBonds()
                        self.drawSelection()
                        if self.parent.cell:
                                self.drawCell()

        def drawAtoms(self):
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
                for atom in self.parent.atoms:
                        #load model matrix with coordinates
                        self.mMatrix.setToIdentity()
                        #move atoms to coord and recognize offset
                        self.mMatrix.translate(atom[0][0],atom[0][1],atom[0][2])
                        #scale atoms depending on species
                        self.mMatrix.scale(atom[2])
                        #bind transformation matrices
                        self.sphereShader.setUniformValue('mvpMatrix',self.proj*self.vMatrix*self.mMatrix)
                        self.sphereShader.setUniformValue('vMatrix',self.vMatrix)
                        self.sphereShader.setUniformValue('mMatrix',self.mMatrix)
                        #create light source:
                        self.sphereShader.setUniformValue('LightPosition_cameraspace',QVector3D(10,10,10))
                        #color vertices
                        self.sphereShader.setUniformValue('MaterialDiffuseColor',atom[1])
                        #glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
                        #draw
                        glDrawArrays(GL_TRIANGLES,0,len(self.atom_modelspace))
                #reset
                self.sphereShader.disableAttributeArray('vertex_modelspace')
                self.sphereShader.disableAttributeArray('normals_modelspace')
                self.sphereShader.release()

        def drawAtomsSelect(self):
                #bind shaders:
                self.selectShader.bind()

                #send vertices
                self.selectShader.setAttributeArray('vertex_modelspace',self.atom_modelspace)
                self.selectShader.enableAttributeArray('vertex_modelspace')

                #render loop
                for i in range(len(self.parent.atoms)):
                        #color from atom number:
                        r = (i & 0x000000FF) >> 0
                        g = (i & 0x0000FF00) >> 8
                        b = (i & 0x00FF0000) >> 16
                        #load model matrix with coordinates
                        self.mMatrix.setToIdentity()
                        atom = self.parent.atoms[i]
                        #move atoms to coord and recognize offset
                        self.mMatrix.translate(atom[0][0],atom[0][1],atom[0][2])
                        #scale atoms depending on species
                        self.mMatrix.scale(atom[2])
                        #bind transformation matrix
                        self.selectShader.setUniformValue('mvpMatrix',self.proj*self.vMatrix*self.mMatrix)
                        #color vertices
                        self.selectShader.setUniformValue('MaterialDiffuseColor',QColor(r,g,b,255))
                        #draw
                        glDrawArrays(GL_TRIANGLES,0,len(self.atom_modelspace))
                #reset
                self.selectShader.disableAttributeArray('vertex_modelspace')
                self.selectShader.disableAttributeArray('normals_modelspace')
                self.selectShader.release()

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
                        self.sphereShader.setUniformValue('vMatrix',self.vMatrix)
                        self.sphereShader.setUniformValue('mMatrix',self.mMatrix)
                        #create light source:
                        self.sphereShader.setUniformValue('LightPosition_cameraspace',QVector3D(10,10,10))
                        #color vertices
                        self.sphereShader.setUniformValue('MaterialDiffuseColor',QColor(120,120,150,200))
                        #glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
                        #draw
                        glDrawArrays(GL_TRIANGLES,0,len(self.atom_modelspace))
                #reset
                self.sphereShader.disableAttributeArray('vertex_modelspace')
                self.sphereShader.disableAttributeArray('normals_modelspace')
                self.sphereShader.release()

        def drawBonds(self):
                self.bondShader.bind()

                #send vertices
                self.bondShader.setAttributeArray('vertex_modelspace',self.bond_modelspace)
                self.bondShader.enableAttributeArray('vertex_modelspace')
                #send normals for lighting
                self.bondShader.setAttributeArray('normals_modelspace',self.bo_normals_modelspace)
                self.bondShader.enableAttributeArray('normals_modelspace')

                #render loop
                for bond in self.parent.bonds:
                        #load model Matrix with coordinates relative to offset
                        self.mMatrix.setToIdentity()
                        self.mMatrix.translate(bond[0][0],bond[0][1],bond[0][2])
                        #rotate to fit bond vector
                        self.mMatrix.rotate(bond[1],bond[2][0],bond[2][1],bond[2][2])
                        # bind transformation matrices
                        self.bondShader.setUniformValue('mvpMatrix',self.proj*self.vMatrix*self.mMatrix)
                        self.bondShader.setUniformValue('vMatrix',self.vMatrix)
                        self.bondShader.setUniformValue('mMatrix',self.mMatrix)
                        #light source:
                        self.bondShader.setUniformValue('LightPosition_cameraspace',QVector3D(10,10,10))
                        #color vertices
                        self.bondShader.setUniformValue('Side1DiffuseColor',bond[3])
                        self.bondShader.setUniformValue('Side2DiffuseColor',bond[4])
                        #draw
                        glDrawArrays(GL_TRIANGLES,0,len(self.bond_modelspace))
                #reset
                self.bondShader.disableAttributeArray('vertex_modelspace')
                self.bondShader.disableAttributeArray('normals_modelspace')
                self.bondShader.release()

        def drawCell(self):
                #bind shaders:
                self.lineShader.bind()
                self.lineShader.setAttributeArray('vertex_modelspace',self.parent.cell_modelspace)
                self.lineShader.enableAttributeArray('vertex_modelspace')
                self.lineShader.setUniformValue('color',QColor(0,0,0))
                for i in self.parent.cells:
                        self.mMatrix.setToIdentity()
                        #move viewpoint to center
                        self.mMatrix.translate(i[0],i[1],i[2])
                        self.lineShader.setUniformValue('mvpMatrix',self.proj*self.vMatrix*self.mMatrix)
                        glDrawArrays(GL_LINES,0,len(self.parent.cell_modelspace))
                self.lineShader.disableAttributeArray('vertex_modelspace')
                self.lineShader.release()

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
                        self.rotMat.setToIdentity()
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
                        self.rotMat = tmp*self.rotMat
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
                        self.rotMat = tmp*self.rotMat
                        self.updateGL()
                if e.key() == Qt.Key_Down:
                        tmp = QMatrix4x4()
                        tmp.rotate(10,1,0,0)
                        self.rotMat = tmp*self.rotMat
                        self.updateGL()
                if e.key() == Qt.Key_Left:
                        tmp = QMatrix4x4()
                        tmp.rotate(-10,0,1,0)
                        self.rotMat = tmp*self.rotMat
                        self.updateGL()
                if e.key() == Qt.Key_Right:
                        tmp = QMatrix4x4()
                        tmp.rotate(10,0,1,0)
                        self.rotMat = tmp*self.rotMat
                        self.updateGL()
