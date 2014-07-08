#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

from os.path import expanduser
from os import getcwd
from functools import partial

from PyQt4.QtGui import *

class CoordTB(QMainWindow):

        def __init__(self,controller):
                super(CoordTB,self).__init__()
                self.controller = controller
                self.initApp()

        def initApp(self):

                #Create Menu
                self.initMenu()

                #Create main widget
                mv = MainView(self.controller)
                self.setCentralWidget(mv)

                # Set Title and run:
                self.setWindowTitle('CoordToolBox')
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
                self.loadAction = QAction('Load',self)
                self.loadAction.setShortcut('Ctrl+O')
                self.loadAction.triggered.connect(self.loadHandler)
                self.saveAction = QAction('Save',self)
                self.saveAction.setShortcut('Ctrl+S')
                self.saveAction.triggered.connect(self.saveHandler)
                self.exitAction = QAction('Exit',self)
                self.exitAction.setShortcut('Ctrl+Q')
                self.exitAction.triggered.connect(qApp.quit)

        def loadHandler(self):
                fname = QFileDialog.getOpenFileName(self,'Open File',getcwd())
                if not fname: return
                ftype = QInputDialog.getItem(self,'Choose file type','File type:',self.controller.indict.keys(),0,False)
                ftype = str(ftype[0])
                self.controller.readFile(ftype,fname)
                self.centralWidget().updateView()

        def saveHandler(self):
                msgBox = QMessageBox()
                msgBox.setText('You had saved your file.\nIf this Program was able to.')
                msgBox.exec_()

class MainView(QWidget):

        def __init__(self,controller):
                super(MainView,self).__init__()
                self.controller = controller
                self.initUI()
                self.tabs = ''

        def initUI(self):

                #Molecules:
                self.mol = EditArea(self.controller)

                self.mlist = QListWidget()
                self.mlist.currentRowChanged.connect(self.mol.setCurrentIndex)

                #PWParameter stuff:
                self.pw = EditArea(self.controller)
                
                self.pwlist = QListWidget()
                self.pwlist.currentRowChanged.connect(self.pw.setCurrentIndex)

                #Layout
                self.vbox = QVBoxLayout()
                self.vbox.addWidget(QLabel('Loaded Molecules:'))
                self.vbox.addWidget(self.mlist)
                self.vbox.addWidget(QLabel('PW Parameter sets:'))
                self.vbox.addWidget(self.pwlist)
                self.hbox = QHBoxLayout()
                self.hbox.addLayout(self.vbox)
                self.hbox.addWidget(self.mol,1)
                self.hbox.addWidget(self.pw,1)

                self.setLayout(self.hbox)

        def updateView(self):
                count = self.mol.count()
                for i in range(count,self.controller.get_nmol()):
                        self.mol.addWidget(MolTab(self.controller.get_mol(i)))#,'Mol '+str(i+1))
                        self.mlist.addItem("Mol "+str(i+1))
                count = self.pw.count()
                for i in range(count,self.controller.get_npw()):
                        self.pw.addWidget(PWTab(self.controller.get_pw(i)))
                        self.pwlist.addItem("Param "+str(i+1))

class EditArea(QStackedWidget):

        def __init__(self,controller):
                super(EditArea,self).__init__()
                self.controller = controller
                self.setMinimumSize(440,350)

class MolTab(QWidget):

        def __init__(self,mol):
                super(MolTab,self).__init__()
                #self.controller = controller
                self.mol = mol
                self.initTab()

        def initTab(self):

                # coord fmt dropdown selector
                self.fmt = QComboBox()
                self.fmt.setToolTip('Select format of coordinates')
                self.fmt.addItem('angstrom')
                self.fmt.addItem('bohr')
                self.fmt.addItem('crystal')
                self.fmt.addItem('alat')
                self.fmt.currentIndexChanged.connect(self.fillTab)

                # show celldm
                self.cellDm = QLineEdit()
                self.cellDm.editingFinished.connect(self.updateCellDm)

                # layout1
                self.hbox = QHBoxLayout()
                self.hbox.addWidget(self.fmt)
                self.hbox.addStretch(1)
                self.hbox.addWidget(QLabel('Cell dimension:'))
                self.hbox.addWidget(self.cellDm)

                # show coordinates in table
                self.table = QTableWidget()
                self.table.setColumnCount(4)
                self.table.setHorizontalHeaderLabels(['Type','x','y','z'])
                self.table.itemChanged.connect(self.cellHandler)

                # show coordinates in table
                self.vtable = QTableWidget()
                self.vtable.setColumnCount(3)
                self.vtable.setRowCount(3)
                self.vtable.setFixedHeight(120)
                self.vtable.setHorizontalHeaderLabels(['x','y','z'])
                self.vtable.itemChanged.connect(self.vecHandler)

                # set Layout for Tab
                self.vbox = QVBoxLayout()
                self.vbox.addLayout(self.hbox)
                self.vbox.addWidget(QLabel('Coordinates:'))
                self.vbox.addWidget(self.table)
                self.vbox.addWidget(QLabel('Cell vectors:'))
                self.vbox.addWidget(self.vtable)

                self.setLayout(self.vbox)
                self.resize(self.sizeHint())

                # initialize content
                self.cellDm.setText(str(self.mol.get_celldm()))
                self.table.setRowCount(self.mol.get_nat())
                self.fmt.setCurrentIndex(0)

                # fill tab with coordinates in Ångström
                self.fillTab()

        def fillTab(self):
                #prevent handling of cell changes during fill
                self.tabledisable = True
                #fill atom table
                for i in range(self.mol.get_nat()):
                        name = QTableWidgetItem(self.mol.get_atom(i).get_name())
                        self.table.setItem(i,0,name)
                        coord = self.mol.get_atom(i).get_coord(self.fmt.currentText())
                        for j in range(3):
                                self.table.setItem(i,j+1,QTableWidgetItem(str(coord[j])))
                #fill cell vec list
                vec = self.mol.get_vec()
                for i in range(3):
                        for j in range(3):
                                self.vtable.setItem(i,j,QTableWidgetItem(str(vec[i,j])))
                #reenable handling
                self.tabledisable = False

        def updateCellDm(self):
                self.mol.set_celldm(float(self.cellDm.text()))
                self.fillTab()

        def cellHandler(self):
                if self.tabledisable: return
                atom = self.table.currentRow()
                if self.table.currentColumn() == 0:
                        self.mol.get_atom(atom).set_name(self.table.item(atom,0).text())
                else:
                        coord = [0,0,0]
                        for j in range(3):
                                coord[j]=float(self.table.item(atom,j+1).text())
                        self.mol.get_atom(atom).set_coord(self.fmt.currentText(),coord)

        def vecHandler(self):
                if self.tabledisable: return
                vec=[[0,0,0],[0,0,0],[0,0,0]]
                for i in range(3):
                        for j in range(3):
                                vec[i][j]=float(self.vtable.item(i,j).text())
                self.mol.set_vec(vec)
                self.fillTab()

class PWTab(QTreeWidget):

        def __init__(self,pw):
                super(PWTab,self).__init__()
                self.pw = pw
                self.setColumnCount(2)
                self.setHeaderLabels(['Parameter','Value'])
                self.makeTree()
                
        def makeTree(self):
                #container for items
                items=[[],[]]
                #mandatory namelists
                for i in ['&control','&system','&electrons']:
                        items[0].append(QTreeWidgetItem(self))
                        items[0][-1].setText(0,i)

                #show optional namelists only if existing
                if '&ions' in self.pw:
                        items[0].append(QTreeWidgetItem(self))
                        items[0][-1].setText(0,'&ions')
                if '&cell' in self.pw:
                        items[0].append(QTreeWidgetItem(self))
                        items[0][-1].setText(0,'&cell')

                #show child entries
                for i in items[0]:
                        for j in self.pw[str(i.text(0))].items():
                                items[1].append(QTreeWidgetItem(i))
                                items[1][-1].setText(0,j[0])
                                items[1][-1].setText(1,j[1])
