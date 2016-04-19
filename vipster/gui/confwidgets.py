# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt
from copy import copy

from vipster.settings import pse,config,saveConfig,default

class ConfWidget(QWidget):
    def __init__(self,parent,data):
        super(ConfWidget,self).__init__()
        self.parent=parent
        self.updateDisable=False
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
        self.updateDisable=True
        self.data = data
        for i in range(len(self.names)):
            if self.types[i] == bool:
                self.widgets[i].setCheckState(2*self.data[self.names[i]])
            elif self.types[i] == list:
                self.widgets[i].setStyleSheet('background-color: #{:02X}{:02X}{:02X}'.format(
                    *[int(255*j) for j in self.data[self.names[i]]]))
            else:
                self.widgets[i].setText(str(self.data[self.names[i]]))
        self.updateDisable=False

    def changeHandler(self,change=None):
        if self.updateDisable:return
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
        self.parent.updateMol()

class Settings(QWidget):
    def __init__(self,parent):
        super(Settings,self).__init__()
        settings = ConfWidget(parent,config)
        settings.setData(config)
        saveBut = QPushButton('Save settings')
        saveBut.clicked.connect(saveConfig)
        loadBut = QPushButton('Load default settings')
        def loadDefault():
            for i in config:
                config[i] = default['General'][i]
            settings.setData(config)
            parent.updateMol()
        loadBut.clicked.connect(loadDefault)
        scroll = QScrollArea()
        scroll.setFrameStyle(0)
        scroll.setWidget(settings)
        layout = QGridLayout()
        layout.addWidget(scroll,0,0,1,2)
        layout.addWidget(saveBut,1,0)
        layout.addWidget(loadBut,1,1)
        self.setLayout(layout)

class PseGlobal(QWidget):
    def __init__(self,parent):
        super(PseGlobal,self).__init__()
        self.parent = parent
        self.settings = ConfWidget(parent,pse['H'])
        self.atomlist = QListWidget()
        self.atomlist.currentTextChanged.connect(self.fillPse)
        self.atomlist.addItems(list(pse.keys()))
        self.atomlist.setCurrentRow(0)
        loadBut = QPushButton('Load default element')
        def loadDefault():
            el = str(self.atomlist.currentItem().text())
            pse[el] = copy(default['PSE'][el])
            self.fillPse(el)
        loadBut.clicked.connect(loadDefault)
        saveBut = QPushButton('Save global pse')
        saveBut.clicked.connect(saveConfig)
        layout = QGridLayout()
        layout.addWidget(self.atomlist,0,0)
        layout.addWidget(self.settings,0,1)
        layout.addWidget(saveBut,1,0)
        layout.addWidget(loadBut,1,1)
        self.setLayout(layout)

    def setMol(self,mol):
        #refresh, in case the global pse is influenced externally
        self.settings.setData(pse[str(self.atomlist.currentItem().text())])

    def fillPse(self,text):
        self.settings.setData(pse[str(text)])

class PseMol(QWidget):
    def __init__(self,parent):
        super(PseMol,self).__init__()
        self.parent = parent
        self.settings = ConfWidget(parent,pse['H'])
        self.atomlist = QListWidget()
        self.atomlist.currentTextChanged.connect(self.fillPse)
        loadBut = QPushButton('Load global element')
        def loadGlobal():
            if not self.atomlist.currentItem(): return
            el = str(self.atomlist.currentItem().text())
            del self.mol.pse[el]
            self.fillPse(el)
            self.parent.updateMol()
        loadBut.clicked.connect(loadGlobal)
        saveBut = QPushButton('Save global element')
        def saveGlobal():
            if not self.atomlist.currentItem(): return
            el = str(self.atomlist.currentItem().text())
            pse[el] = copy(self.mol.pse[el])
            self.parent.updateMol()
        saveBut.clicked.connect(saveGlobal)
        layout = QGridLayout()
        layout.addWidget(self.atomlist,0,0)
        layout.addWidget(self.settings,0,1)
        layout.addWidget(saveBut,1,0)
        layout.addWidget(loadBut,1,1)
        self.setLayout(layout)

    def setMol(self,mol):
        if hasattr(self,'mol') and self.mol is mol: return
        self.mol = mol
        self.atomlist.clear()
        self.atomlist.addItems(list(mol.pse.keys()))

    def fillPse(self,text):
        if hasattr(self,'mol'):
            self.settings.setData(self.mol.pse[str(text)])
