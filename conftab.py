#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from copy import copy

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

        pbox = QVBoxLayout()
        self.pseEdits = []
        for i in ['Pseudopotential:','Z:','Mass:','Covalent radius:','VdW radius:']:
            pbox.addWidget(QLabel(i))
            self.pseEdits.append(QLineEdit())
            pbox.addWidget(self.pseEdits[-1])
            self.pseEdits[-1].editingFinished.connect(self.elemHandler)
        self.pseEdits[1].setValidator(QIntValidator())
        self.pseEdits[2].setValidator(QDoubleValidator())
        self.pseEdits[3].setValidator(QDoubleValidator())
        self.pseEdits[4].setValidator(QDoubleValidator())
        pbox.addWidget(QLabel('Color:'))
        self.color = QColor()
        self.colBut = QPushButton()
        self.colBut.clicked.connect(self.colHandler)
        pbox.addWidget(self.colBut)
        pbox.addStretch()
        self.pseWidget = QWidget()
        self.pseWidget.setLayout(pbox)

        self.settings = self.parent.controller.config
        self.setWidget = QWidget()
        self.boolLabel = QLabel('Enabled:')
        self.boolWidget = QCheckBox()
        self.boolWidget.stateChanged.connect(self.boolHandler)
        self.valLabel = QLabel('Value:')
        self.valWidget = QLineEdit()
        self.valWidget.editingFinished.connect(self.valHandler)
        sbox=QVBoxLayout()
        sbox.addWidget(self.boolLabel)
        sbox.addWidget(self.boolWidget)
        sbox.addWidget(self.valLabel)
        sbox.addWidget(self.valWidget)
        sbox.addStretch()
        self.setWidget.setLayout(sbox)

        self.fileBut = QPushButton('Save config')
        self.fileBut.clicked.connect(self.saveFile)
        self.setDefBut = QPushButton('Load default setting')
        self.setDefBut.clicked.connect(self.settingFromDefault)
        self.pseDefBut = QPushButton('Load default element')
        self.pseDefBut.clicked.connect(self.pseFromDefault)
        self.pseLoadBut = QPushButton('Load global element')
        self.pseLoadBut.clicked.connect(self.pseFromInstance)
        self.pseSavBut = QPushButton('Save element to global')
        self.pseSavBut.clicked.connect(self.pseToInstance)

        hbox = QHBoxLayout()
        hbox.addWidget(self.list)
        hbox.addWidget(self.pseWidget)
        hbox.addWidget(self.setWidget)
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.fileBut)
        hbox2.addWidget(self.setDefBut)
        hbox2.addWidget(self.pseDefBut)
        hbox2.addWidget(self.pseSavBut)
        hbox2.addWidget(self.pseLoadBut)
        vbox = QVBoxLayout()
        vbox.addWidget(self.selection)
        vbox.addLayout(hbox)
        vbox.addLayout(hbox2)
        self.setLayout(vbox)

    def setMol(self,mol):
        if hasattr(self,'mol') and self.mol is mol: return
        self.mol = mol
        self.fillTab()

    def fillTab(self):
        sel = self.selection.currentText()
        if sel == 'Settings':
            self.pseWidget.hide()
            self.setWidget.show()
            self.fileBut.show()
            self.setDefBut.show()
            self.pseDefBut.hide()
            self.pseLoadBut.hide()
            self.pseSavBut.hide()
            self.list.clear()
            self.list.addItems(list(self.settings.keys()))
            self.list.currentTextChanged.disconnect()
            self.list.currentTextChanged.connect(self.setSetting)
        elif sel == 'PSE':
            self.setWidget.hide()
            self.pseWidget.show()
            self.fileBut.show()
            self.setDefBut.hide()
            self.pseDefBut.show()
            self.pseLoadBut.hide()
            self.pseSavBut.hide()
            self.list.clear()
            self.pse = self.parent.controller.pse
            self.list.addItems(list(self.pse.keys()))
            self.list.currentTextChanged.disconnect()
            self.list.currentTextChanged.connect(self.setElement)
        elif sel == 'PSE (Mol)':
            self.setWidget.hide()
            self.pseWidget.show()
            self.fileBut.hide()
            self.setDefBut.hide()
            self.pseDefBut.hide()
            self.pseLoadBut.show()
            self.pseSavBut.show()
            self.list.clear()
            self.pse = self.mol.pse
            self.list.addItems(list(self.pse.keys()))
            self.list.currentTextChanged.disconnect()
            self.list.currentTextChanged.connect(self.setElement)
        self.list.setCurrentRow(0)

    def setElement(self,e):
        if not e:
            self.element = None
            for i in range(5):
                self.pseEdits[i].setText('')
                self.pseEdits[i].setDisabled(True)
            self.colBut.setStyleSheet('')
            self.colBut.setDisabled(True)
        else:
            self.element=str(e)
            e = self.pse[self.element]
            for i in range(5):
                self.pseEdits[i].setText(str(e[i]))
                self.pseEdits[i].setEnabled(True)
            self.color.setRgbF(*e[5:])
            self.colBut.setEnabled(True)
            self.colBut.setStyleSheet('background-color: #{:06X}'.format(
                    self.color.rgb()&0xFFFFFF))

    def setSetting(self,s):
        if not s:
            self.setting = None
            return
        else:
            self.setting = str(s)
            val = self.settings[self.setting]
            if type(val) is bool:
                self.valLabel.setDisabled(True)
                self.valWidget.setDisabled(True)
                self.valWidget.setText('')
                self.boolLabel.setEnabled(True)
                self.boolWidget.setEnabled(True)
                if val:
                    self.boolWidget.setChecked(True)
                else:
                    self.boolWidget.setChecked(False)
            elif type(val) is float:
                self.boolLabel.setDisabled(True)
                self.boolWidget.setDisabled(True)
                self.boolWidget.setChecked(False)
                self.valLabel.setEnabled(True)
                self.valWidget.setEnabled(True)
                self.valWidget.setValidator(QDoubleValidator())
                self.valWidget.setText(str(val))
            else:
                self.boolLabel.setDisabled(True)
                self.boolWidget.setDisabled(True)
                self.boolWidget.setChecked(False)
                self.valLabel.setEnabled(True)
                self.valWidget.setEnabled(True)
                self.valWidget.setValidator(None)
                self.valWidget.setText(str(val))

    def elemHandler(self):
        e = self.pse[self.element]
        e[0] = self.pseEdits[0].text()
        e[1] = int(self.pseEdits[1].text())
        e[2] = float(self.pseEdits[2].text())
        e[3] = float(self.pseEdits[3].text())
        e[4] = float(self.pseEdits[4].text())
        self.parent.updateMolStep()

    def colHandler(self):
        col = QColorDialog.getColor(self.color,self,'Select Color',QColorDialog.ShowAlphaChannel)
        if col.isValid():
            e=self.pse[self.element]
            e[5]=col.redF()
            e[6]=col.greenF()
            e[7]=col.blueF()
            e[8]=col.alphaF()
        self.setElement(self.element)
        self.parent.updateMolStep()

    def boolHandler(self,check):
        if check:
            self.settings[self.setting] = True
        else:
            self.settings[self.setting] = False
        self.parent.updateMolStep()

    def valHandler(self):
        if self.valWidget.validator():
            self.settings[self.setting] = float(self.valWidget.text())
        else:
            self.settings[self.setting] = str(self.valWidget.text())
        self.parent.updateMolStep()

    def saveFile(self):
        self.parent.controller.saveConfig()

    def settingFromDefault(self):
        self.settings[self.setting] = copy(self.parent.controller.default['General'][self.setting])
        self.setSetting(self.setting)

    def pseFromDefault(self):
        self.pse[self.element] = copy(self.parent.controller.default['PSE'][self.element])
        self.valHandler()
        self.setElement(self.element)

    def pseFromInstance(self):
        if self.element:
            del self.pse[self.element]
            self.setElement(self.element)
            self.parent.updateMolStep()

    def pseToInstance(self):
        self.parent.controller.pse[self.element] = copy(self.pse[self.element])
