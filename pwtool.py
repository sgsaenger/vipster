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

        def initUI(self):

                #add main Views
                iv = Input(self.controller)
                ov = Output(self.controller)
                mod = Modify()
                tab = EditArea(self.controller)

                #Layout
                vbox = QVBoxLayout()
                vbox.addWidget(iv)
                vbox.addWidget(mod)
                vbox.addWidget(ov)
                vbox.addStretch()
                hbox = QHBoxLayout()
                hbox.addLayout(vbox)
                hbox.addWidget(tab,1)

                self.setLayout(hbox)


class Input(QGroupBox):

        def __init__(self,controller):
                super(Input,self).__init__()
                self.controller = controller
                self.initBox()

        def initBox(self):
                #Open button
                inpOpen = QPushButton('Open',self)
                inpOpen.setToolTip('Search for file to open')
                inpOpen.resize(inpOpen.sizeHint())
                inpOpen.clicked.connect(self.openDiag)

                #text field with path to file
                self.inpPath = QLineEdit()

                #Filetype selector
                self.inpSel = QComboBox(self)
                self.inpSel.setToolTip('Select format of input file')
                self.inpSel.addItem("PWScf Input")
                self.inpSel.addItem("PWScf Output")
                self.inpSel.addItem("xyz")

                #Load button
                inpLoad = QPushButton('Load',self)
                inpLoad.setToolTip('Load input file of specified format')
                inpLoad.resize(inpLoad.sizeHint())
                inpLoad.clicked.connect(self.loadHandler)


                #Layout. Good.
                vbox = QVBoxLayout()
                hbox = QHBoxLayout()
                hbox.addWidget(self.inpSel)
                hbox.addStretch()
                hbox.addWidget(inpOpen)
                hbox.addWidget(inpLoad)
                vbox.addWidget(self.inpPath)
                vbox.addLayout(hbox)


                #Finalize inpView
                self.setLayout(vbox)
                self.setTitle('Input')

        def openDiag(self):
                fname = QFileDialog.getOpenFileName(self, 'Open file',getcwd())
                self.inpPath.setText(fname)
        def loadHandler(self):
                self.controller.readFile(str(self.inpSel.currentText()),str(self.inpPath.text()))


class Output(QGroupBox):

        def __init__(self,controller):
                super(Output,self).__init__()
                self.controller = controller
                self.initBox()

        def initBox(self):

                #Output selector drop down and save button
                outSave = QPushButton('Save',self)
                outSave.setToolTip('Save in specified format')
                outSave.resize(outSave.sizeHint())

                #text field with path to file
                self.outPath = QLineEdit()

                outSel = QComboBox(self)
                outSel.setToolTip('Select format of output file')
                outSel.addItem("PWScf Input")
                outSel.addItem("xyz")

                #Open button
                outOpen = QPushButton('Open',self)
                outOpen.setToolTip('Search for file to open')
                outOpen.resize(outOpen.sizeHint())
                outOpen.clicked.connect(self.saveDiag)

                vbox = QVBoxLayout()
                hbox = QHBoxLayout()
                hbox.addWidget(outSel)
                hbox.addStretch()
                hbox.addWidget(outOpen)
                hbox.addWidget(outSave)
                vbox.addWidget(self.outPath)
                vbox.addLayout(hbox)

                self.setTitle('Output')
                self.setLayout(vbox)

        def saveDiag(self):
                fname = QFileDialog.getSaveFileName(self, 'Open file',expanduser("~"))
                self.outPath.setText(fname)

        def saveHandler(self):
                f = open(self.outPath.text(),'w')
                with f:
                        f.write('here be dragons')



class EditArea(QTabWidget):

        def __init__(self,controller):
                super(EditArea,self).__init__()
                self.controller = controller
                self.initArea()

        def initArea(self):

                #Table with coordinates, different formats (see xcrysden?)
                coordTab = QWidget()

                #CellDM and format of coordinates
                coordFmt = QComboBox()
                coordFmt.setToolTip('Select format of coordinates')
                coordFmt.addItem('Angstrom')
                coordFmt.addItem('Bohr')
                coordFmt.addItem('Crystal')
                coordFmt.addItem('Alat')

                cellDm = QLineEdit()

                hbox = QHBoxLayout()
                hbox.addWidget(coordFmt)
                hbox.addStretch(1)
                hbox.addWidget(QLabel('Cell dimension:'))
                hbox.addWidget(cellDm)

                #Table
                coordTable = QTableWidget()
                coordTable.setColumnCount(5)
                coordTable.setHorizontalHeaderLabels(['#','Type','x','y','z'])

                #set Layout for Tab
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(coordTable)

                coordTab.setLayout(vbox)

                #Table or similar with calculation parameters
                Calc = QTableWidget()

                self.addTab(coordTab,'Coordinates')
                self.addTab(Calc,'Calc. param.')

class Modify(QGroupBox):

        def __init__(self):
                super(Modify,self).__init__()
                self.initMod()

        def initMod(self):

                self.setCheckable(True)
                self.setChecked(False)
                self.setTitle('Modify')

def main():
        app = QApplication(sys.argv)
        o = ''
        main = CoordTB(o)
        sys.exit(app.exec_())

