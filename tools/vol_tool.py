#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

####################################
# Volume-PP
####################################

class Volume(QWidget):

    def __init__(self,parent):
        super(Volume,self).__init__()

        def volSHandler():
            if hasattr(self,'vol'):
                curVal = self.volSel.value()/1000. * (self.vol.max()-self.vol.min()) + self.vol.min()
                self.volCur.setText(str(curVal))
                parent.visual.setSurf(curVal)
            else:
                self.volCur.setText('0')

        def volBHandler():
            volSHandler()
            parent.visual.toggleSurf()

        self.volMin = QLabel('0')
        self.volSel = QSlider()
        self.volSel.setDisabled(True)
        self.volSel.setOrientation(1)
        self.volSel.setMinimum(0)
        self.volSel.setMaximum(0)
        self.volSel.setTickPosition(self.volSel.TicksBelow)
        self.volSel.setSingleStep(1)
        self.volSel.setTickInterval(100)
        self.volSel.valueChanged.connect(volSHandler)
        self.volCur = QLabel('0')
        self.volBut = QPushButton('Show/Hide')
        self.volBut.clicked.connect(volBHandler)
        self.volBut.setDisabled(True)
        self.volMax = QLabel('0')
        vbox = QVBoxLayout()
        hbox = QHBoxLayout()
        hbox.addWidget(QLabel('Range:'))
        hbox.addWidget(self.volMin)
        hbox.addWidget(QLabel(' to '))
        hbox.addWidget(self.volMax)
        vbox.addLayout(hbox)
        vbox.addWidget(self.volSel)
        vbox.addWidget(self.volCur)
        vbox.addWidget(self.volBut)
        vbox.addStretch()
        self.setLayout(vbox)

    def setMol(self,mol):
        if hasattr(mol,'_vol'):
            self.vol = mol.get_vol()
            self.volMin.setText(str(self.vol.min()))
            self.volSel.setMaximum(1000)
            self.volMax.setText(str(self.vol.max()))
            self.volSel.setEnabled(True)
            self.volBut.setEnabled(True)
        else:
            self.volMin.setText('0')
            self.volSel.setMaximum(0)
            self.volMax.setText('0')
            self.volSel.setDisabled(True)
            self.volBut.setDisabled(True)

    def setSel(self,sel):
        return
