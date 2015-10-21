# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt
from copy import copy

class ConfTab(QWidget):
    def __init__(self,parent):
        super(ConfTab,self).__init__()
        self.parent = parent

        self.selection = QComboBox()
        for i in ['PSE (Mol)','PSE','Settings']:
            self.selection.addItem(i)
        self.selection.currentIndexChanged.connect(self.fillTab)

        self.setWidget = ConfTab.ConfWidget(parent,parent.controller.config)
        self.setWidget.setData(parent.controller.config)
        self.setWidget.hide()

        self.pseList = QListWidget()
        self.pseList.currentTextChanged.connect(self.fillPse)
        self.pseArea = ConfTab.ConfWidget(parent,parent.controller.pse['X'])
        psebox = QHBoxLayout()
        psebox.addWidget(self.pseList)
        psebox.addWidget(self.pseArea)
        self.pseWidget = QWidget()
        self.pseWidget.setLayout(psebox)

        self.setDefaultBut = QPushButton('Load default settings')
        self.setDefaultBut.clicked.connect(self.configFromDefault)
        self.setFileBut = QPushButton('Save config')
        self.setFileBut.clicked.connect(self.saveFile)
        self.pseFileBut = QPushButton('Save config')
        self.pseFileBut.clicked.connect(self.saveFile)
        self.pseDefaultBut = QPushButton('Load default element')
        self.pseDefaultBut.clicked.connect(self.elementFromDefault)
        self.pseLoadBut = QPushButton('Load global element')
        self.pseLoadBut.clicked.connect(self.elementFromGlobal)
        self.pseSaveBut = QPushButton('Save global element')
        self.pseSaveBut.clicked.connect(self.elementToGlobal)
        buttongrid = QGridLayout()
        buttongrid.addWidget(self.setDefaultBut,0,0)
        buttongrid.addWidget(self.setFileBut,0,1)
        buttongrid.addWidget(self.pseDefaultBut,1,0)
        buttongrid.addWidget(self.pseFileBut,1,1)
        buttongrid.addWidget(self.pseLoadBut,2,0)
        buttongrid.addWidget(self.pseSaveBut,2,1)

        vbox = QVBoxLayout()
        vbox.addWidget(self.selection)
        vbox.addWidget(self.setWidget)
        vbox.addWidget(self.pseWidget)
        vbox.addLayout(buttongrid)
        self.setLayout(vbox)

    def saveFile(self):
        self.parent.controller.saveConfig()

    def configFromDefault(self):
        self.parent.controller.config = copy(self.parent.controller.default['General'])
        self.setWidget.setData(self.parent.controller.config)
        self.parent.updateMolStep()

    def elementFromDefault(self):
        element = str(self.pseList.currentItem().text())
        self.parent.controller.pse[element] = copy(self.parent.controller.default['PSE'][element])
        self.fillPse(element)

    def elementFromGlobal(self):
        if not self.pseList.currentItem(): return
        element = str(self.pseList.currentItem().text())
        del self.mol.pse[element]
        self.fillPse(element)
        self.parent.updateMolStep()

    def elementToGlobal(self):
        if not self.pseList.currentItem(): return
        element = str(self.pseList.currentItem().text())
        self.parent.controller.pse[element] = copy(self.pse[element])

    def fillTab(self):
        sel = self.selection.currentText()
        if sel == 'Settings':
            self.setWidget.show()
            self.pseWidget.hide()
            self.setDefaultBut.show()
            self.setFileBut.show()
            self.pseDefaultBut.hide()
            self.pseFileBut.hide()
            self.pseLoadBut.hide()
            self.pseSaveBut.hide()
        elif sel == 'PSE':
            self.pse = self.parent.controller.pse
            self.setWidget.hide()
            self.pseWidget.show()
            self.setDefaultBut.hide()
            self.setFileBut.hide()
            self.pseDefaultBut.show()
            self.pseFileBut.show()
            self.pseLoadBut.hide()
            self.pseSaveBut.hide()
            self.pseList.clear()
            self.pseList.addItems(list(self.pse.keys()))
        elif sel == 'PSE (Mol)':
            self.pse = self.mol.pse
            self.setWidget.hide()
            self.pseWidget.show()
            self.setDefaultBut.hide()
            self.setFileBut.hide()
            self.pseDefaultBut.hide()
            self.pseFileBut.hide()
            self.pseLoadBut.show()
            self.pseSaveBut.show()
            self.pseList.clear()
            self.pseList.addItems(list(self.pse.keys()))
        self.pseList.setCurrentRow(0)

    def fillPse(self,text):
        if text:
            self.pseArea.setData(self.pse[str(text)])

    def setMol(self,mol):
        if hasattr(self,'mol') and self.mol is mol: return
        self.mol = mol
        self.fillTab()

    class ConfWidget(QWidget):
        def __init__(self,parent,data):
            super(ConfTab.ConfWidget,self).__init__()
            self.parent=parent
            vbox = QVBoxLayout()
            self.names=[]
            self.widgets=[]
            self.types=[]
            self.color = QColor()
            for i in data.items():
                self.names.append(i[0])
                vbox.addWidget(QLabel(i[0]+':'))
                self.types.append(type(i[1]))
                if isinstance(i[1],bool):
                    self.widgets.append(QCheckBox())
                    self.widgets[-1].stateChanged.connect(self.changeHandler)
                elif isinstance(i[1],int):
                    self.widgets.append(QLineEdit())
                    self.widgets[-1].editingFinished.connect(self.changeHandler)
                    self.widgets[-1].setValidator(QIntValidator())
                elif isinstance(i[1],float):
                    self.widgets.append(QLineEdit())
                    self.widgets[-1].setValidator(QDoubleValidator())
                    self.widgets[-1].editingFinished.connect(self.changeHandler)
                elif isinstance(i[1],list):
                    self.widgets.append(QPushButton())
                    self.widgets[-1].clicked.connect(self.changeHandler)
                else:
                    self.widgets.append(QLineEdit())
                    self.widgets[-1].editingFinished.connect(self.changeHandler)
                vbox.addWidget(self.widgets[-1])
            vbox.addStretch()
            self.setLayout(vbox)

        def setData(self,data):
            self.data = data
            for i in range(len(self.names)):
                if self.types[i] == bool:
                    self.widgets[i].setCheckState(2*self.data[self.names[i]])
                elif self.types[i] == list:
                    self.widgets[i].setStyleSheet('background-color: #{:02X}{:02X}{:02X}'.format(
                        *[int(255*j) for j in self.data[self.names[i]]]))
                else:
                    self.widgets[i].setText(str(self.data[self.names[i]]))

        def changeHandler(self,change=None):
            index = self.widgets.index(self.sender())
            iwidg = self.widgets[index]
            iname = self.names[index]
            itype = self.types[index]
            if itype == bool:
                self.data[iname] = itype(change)
            elif itype == list:
                c = QColorDialog.getColor(QColor(*[int(255*i) for i in self.data[iname]]),
                    self,'Select Color',QColorDialog.ShowAlphaChannel)
                if c.isValid():
                    self.data[iname] = [c.redF(),c.greenF(),c.blueF(),c.alphaF()]
            else:
                self.data[iname] = itype(iwidg.text())
            self.setData(self.data)
            self.parent.updateMolStep()
