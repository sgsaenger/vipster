#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

from ptb_mol import *
from os.path import expanduser

from PyQt4.QtGui import *

class CoordTB(QMainWindow):

        def __init__(self):
                super(CoordTB,self).__init__()
                self.initApp()
                self.controller = TBController()

        def initApp(self):

                #Just define central widget and window settings for now
                mv = MainView()
                self.setCentralWidget(mv)
                #self.setGeometry(300,300,250,150)
                self.setWindowTitle('CoordToolBox')
                self.show()

        def parsePwi(self):
                msgBox = QMessageBox()
                msgBox.setText('Loading PWScf Input file')
                msgBox.exec_()

        def parsePwo(self):
                msgBox = QMessageBox()
                msgBox.setText("Loading PWScf Output file")
                msgBox.exec_()

        def parseXyz(self):
                msgBox = QMessageBox()
                msgBox.setText("Loading xyz file")
                msgBox.exec_()

        def parseFile(self,fmt,data):
                options = {0 : self.parsePwi,
                           1 : self.parsePwo,
                           2 : self.parseXyz,
                          }
                options[fmt]()

        def exportPwi(self):
                msgBox = QMessageBox()
                msgBox.setText('Saving PWScf Input file')
                msgBox.exec_()

        def exportPwo(self):
                msgBox = QMessageBox()
                msgBox.setText("Saving PWScf Output file")
                msgBox.exec_()

        def exportXyz(self):
                msgBox = QMessageBox()
                msgBox.setText("Saving xyz file")
                msgBox.exec_()

        def exportFile(self,fmt,data):
                options = {0 : self.exportPwi,
                           1 : self.exportXyz,
                          }
                options[fmt]()


class MainView(QWidget):

        def __init__(self):
                super(MainView,self).__init__()
                self.initUI()

        def initUI(self):

                #add main Views
                iv = Input()
                ov = Output()
                mod = Modify()
                tab = EditArea()

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

        def __init__(self):
                super(Input,self).__init__()
                self.initBox()

        def initBox(self):
                #Open button
                inpOpen = QPushButton('Open',self)
                inpOpen.setToolTip('Search for file to open')
                inpOpen.resize(inpOpen.sizeHint())
                inpOpen.clicked.connect(self.openDiag)

                #text field with path to file
                self.inpPath = QLineEdit()

                #Load button
                inpLoad = QPushButton('Load',self)
                inpLoad.setToolTip('Load input file of specified format')
                inpLoad.resize(inpLoad.sizeHint())
                inpLoad.clicked.connect(self.loadHandler)

                #Filetype selector
                self.inpSel = QComboBox(self)
                self.inpSel.setToolTip('Select format of input file')
                self.inpSel.addItem("PWScf Input")
                self.inpSel.addItem("PWScf Output")
                self.inpSel.addItem("xyz")


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
                fname = QFileDialog.getOpenFileName(self, 'Open file',expanduser("~"))
                self.inpPath.setText(fname)
        def loadHandler(self):
                f = open(self.inpPath.text(),'r')
                with f:
                        data = f.read()
                self.window().parseFile(self.inpSel.currentIndex(),data)


class Output(QGroupBox):

        def __init__(self):
                super(Output,self).__init__()
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

        def __init__(self):
                super(EditArea,self).__init__()
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
                hbox.addWidget(QLabel('Celldimension:'))
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
        main = CoordTB()
        sys.exit(app.exec_())

if __name__ == '__main__':
        main()
