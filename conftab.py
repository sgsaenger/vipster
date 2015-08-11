#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

class ConfTab(QWidget):
    def __init__(self,parent):
        super(ConfTab,self).__init__()
        self.parent = parent
        self.selection = QComboBox()
        for i in ['PSE (Mol)','PSE','Settings']:
            self.selection.addItem(i)
        self.selection.currentIndexChanged.connect(self.fillTab)
        self.list = QListWidget()
        self.list.currentTextChanged.connect(self.setElement)
        self.pseWidget = self.makePseWidget()
        self.settings = self.parent.controller._config
        self.setWidget = self.makeSetWidget()
        hbox = QHBoxLayout()
        hbox.addWidget(self.list)
        hbox.addWidget(self.pseWidget)
        hbox.addWidget(self.setWidget)
        vbox = QVBoxLayout()
        vbox.addWidget(self.selection)
        vbox.addLayout(hbox)
        self.setLayout(vbox)

    def fillTab(self):
        sel = self.selection.currentText()
        if sel == 'Settings':
            self.pseWidget.hide()
            self.setWidget.show()
            self.list.clear()
            self.list.addItems(self.settings.keys())
            self.list.currentTextChanged.disconnect()
            self.list.currentTextChanged.connect(self.setSetting)
        elif sel == 'PSE':
            self.setWidget.hide()
            self.pseWidget.show()
            self.list.clear()
            self.pse = self.parent.controller.pse
            self.list.addItems(self.pse.keys())
            self.list.currentTextChanged.disconnect()
            self.list.currentTextChanged.connect(self.setElement)
        elif sel == 'PSE (Mol)':
            self.setWidget.hide()
            self.pseWidget.show()
            self.list.clear()
            self.pse = self.mol.pse
            self.list.addItems(self.pse.keys())
            self.list.currentTextChanged.disconnect()
            self.list.currentTextChanged.connect(self.setElement)

    def makePseWidget(self):
        vbox = QVBoxLayout()
        self.pseEdits = []
        for i in ['Pseudopotential:','Z:','Mass:','Covalent radius:','VdW radius:']:
            vbox.addWidget(QLabel(i))
            self.pseEdits.append(QLineEdit())
            vbox.addWidget(self.pseEdits[-1])
            self.pseEdits[-1].editingFinished.connect(self.saveElement)
        self.pseEdits[1].setValidator(QIntValidator())
        self.pseEdits[2].setValidator(QDoubleValidator())
        self.pseEdits[3].setValidator(QDoubleValidator())
        self.pseEdits[4].setValidator(QDoubleValidator())
        vbox.addWidget(QLabel('Color:'))
        self.color = QColor()
        self.colBut = QPushButton()
        self.colBut.clicked.connect(self.colHandler)
        vbox.addWidget(self.colBut)
        vbox.addStretch()
        element = QWidget()
        element.setLayout(vbox)
        return element

    def setElement(self,e):
        if not e: return
        e = self.pse[str(e)]
        for i in range(5):
            self.pseEdits[i].setText(str(e[i]))
        self.color.setRgbF(*e[5:])
        self.colBut.setStyleSheet('background-color: #{:06X}'.format(
                self.color.rgb()&0xFFFFFF))

    def saveElement(self):
        e = self.pse[str(self.list.currentItem().text())]
        e[0] = self.pseEdits[0].text()
        e[1] = int(self.pseEdits[1].text())
        e[2] = float(self.pseEdits[2].text())
        e[3] = float(self.pseEdits[3].text())
        e[4] = float(self.pseEdits[4].text())
        self.parent.updateMolStep()

    def colHandler(self):
        col = QColorDialog.getColor(self.color,self,'Select Color',QColorDialog.ShowAlphaChannel)
        if col.isValid():
            e=self.pse[str(self.list.currentItem().text())]
            e[5]=col.redF()
            e[6]=col.greenF()
            e[7]=col.blueF()
            e[8]=col.alphaF()
        self.setElement(self.list.currentItem().text())
        self.parent.updateMolStep()

    def makeSetWidget(self):
        settings = QWidget()
        self.boolLabel = QLabel('Enabled:')
        self.boolWidget = QCheckBox()
        self.boolWidget.stateChanged.connect(self.boolHandler)
        self.floatLabel = QLabel('Value:')
        self.floatWidget = QLineEdit()
        self.floatWidget.setValidator(QDoubleValidator())
        self.floatWidget.editingFinished.connect(self.valHandler)
        vbox=QVBoxLayout()
        vbox.addWidget(self.boolLabel)
        vbox.addWidget(self.boolWidget)
        vbox.addWidget(self.floatLabel)
        vbox.addWidget(self.floatWidget)
        vbox.addStretch()
        settings.setLayout(vbox)
        return settings

    def setSetting(self,s):
        val = self.settings[str(s)]
        if type(val) is bool:
            self.floatLabel.hide()
            self.floatWidget.hide()
            self.boolLabel.show()
            self.boolWidget.show()
            if val:
                self.boolWidget.setChecked(True)
            else:
                self.boolWidget.setChecked(False)
        elif type(val) is float:
            self.boolLabel.hide()
            self.boolWidget.hide()
            self.floatLabel.show()
            self.floatWidget.show()
            self.floatWidget.setText(str(val))
        pass

    def boolHandler(self,check):
        if check:
            self.settings[str(self.list.currentItem().text())] = True
        else:
            self.settings[str(self.list.currentItem().text())] = False
        self.parent.updateMolStep()

    def valHandler(self):
        self.settings[str(self.list.currentItem().text())] = float(self.floatWidget.text())
        self.parent.updateMolStep()

    def setMol(self,mol):
        if hasattr(self,'mol') and self.mol is mol: return
        self.mol = mol
        self.fillTab()
