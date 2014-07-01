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

                # show celldm
                self.cellDm = QLineEdit()

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

                # set Layout for Tab
                self.vbox = QVBoxLayout()
                self.vbox.addLayout(self.hbox)
                self.vbox.addWidget(self.table)

                self.setLayout(self.vbox)

        def fillTab(self):
                self.cellDm.setText(str(self.mol.get_celldm()))
                self.table.setRowCount(self.mol.get_nat())
                self.fmt.setCurrentIndex(0)
                for i in range(self.mol.get_nat()):
                        name = QTableWidgetItem(self.mol.get_atom(i).get_name())
                        self.table.setItem(i,0,name)
                        coord = self.mol.get_atom(i).get_coord(self.fmt.currentText())
                        for j in range(3):
                                self.table.setItem(i,j+1,QTableWidgetItem(str(coord[j])))

        def updateTab(self):
                self.mol.set_celldm(float(self.cellDm.text()))

class MolArea(QTabWidget):

        def __init__(self,controller):
                super(MolArea,self).__init__()
                self.controller = controller

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
                fname = QFileDialog.getOpenFileName(self,'Open file',getcwd())
                ftype = QInputDialog.getItem(self,'Choose file type','File type:',self.controller.indict.keys(),0,False)
                ftype = str(ftype[0])
                self.controller.readFile(ftype,fname)
                self.parent().updateView()
        def saveHandler(self):
                msgBox = QMessageBox()
                msgBox.setText('You had saved your file.\nIf this Program was able to.')
                msgBox.exec_()
