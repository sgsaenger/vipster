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

                #Just define central widget and window settings for now
                mv = MainView(self.controller)
                self.setCentralWidget(mv)
                self.setWindowTitle('CoordToolBox')
                self.show()

class MainView(QWidget):

        def __init__(self,controller):
                super(MainView,self).__init__()
                self.controller = controller
                self.initUI()
                self.tabs = ''

        def initUI(self):

                #add main Views
                self.io = IOBox(self.controller)
                self.mol = MolArea(self.controller)

                #Layout
                self.vbox = QVBoxLayout()
                self.vbox.addWidget(self.io)
                self.vbox.addStretch()
                self.hbox = QHBoxLayout()
                self.hbox.addLayout(self.vbox)
                self.hbox.addWidget(self.mol,1)

                self.setLayout(self.hbox)

        def updateView(self):
                count = self.mol.count()
                for i in range(count,len(self.controller.mol)):
                        self.mol.addTab(MolTab(self.controller.mol[i]),'Mol '+str(i+1))
                        self.mol.widget(i).fillTab()

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
                self.fmt.addItem(u'Ångström')
                self.fmt.addItem('Bohr')
                self.fmt.addItem('Crystal')
                self.fmt.addItem('Alat')
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
                self.tabledisable = True
                for i in range(self.mol.get_nat()):
                        name = QTableWidgetItem(self.mol.get_atom(i).get_name())
                        self.table.setItem(i,0,name)
                        coord = self.mol.get_atom(i).get_coord(self.fmt.currentText())
                        for j in range(3):
                                self.table.setItem(i,j+1,QTableWidgetItem(str(coord[j])))
                vec = self.mol.get_vec()
                for i in range(3):
                        for j in range(3):
                                self.vtable.setItem(i,j,QTableWidgetItem(str(vec[i,j])))
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

class MolArea(QTabWidget):

        def __init__(self,controller):
                super(MolArea,self).__init__()
                self.controller = controller
                self.setMinimumSize(440,350)

class IOBox(QGroupBox):

        def __init__(self,controller):
                super(IOBox,self).__init__()
                self.controller = controller
                self.initBox()

        def initBox(self):
                #Load button
                load = QPushButton('Load',self)
                load.setToolTip('Load input file')
                load.resize(load.sizeHint())
                load.clicked.connect(self.loadHandler)

                #Output selector drop down and save button
                save = QPushButton('Save',self)
                save.setToolTip('Save in specified format')
                save.resize(save.sizeHint())
                save.clicked.connect(self.saveHandler)

                #Layout. Good.
                hbox = QHBoxLayout()
                hbox.addWidget(load)
                hbox.addWidget(save)

                #Finalize inpView
                self.setLayout(hbox)
                self.setTitle('I/O')

        def loadHandler(self):
                fname = QFileDialog.getOpenFileName(self,'Open File',getcwd())
                ftype = QInputDialog.getItem(self,'Choose file type','File type:',self.controller.indict.keys(),0,False)
                ftype = str(ftype[0])
                self.controller.readFile(ftype,fname)
                self.parent().updateView()
        def saveHandler(self):
                msgBox = QMessageBox()
                msgBox.setText('You had saved your file.\nIf this Program was able to.')
                msgBox.exec_()

        def saveProto(self):
                fname = QFileDialog.getSaveFileName(self,'Save File',getcwd())
