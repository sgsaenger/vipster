#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

####################################
# Crystal and Volume planes
####################################

class Plane(QWidget):

    def __init__(self,parent):
        super(Plane,self).__init__()

        def crysSelect():
            parent.visual.setPlane('c',[ph.value(),pk.value(),pl.value()])

        def crysButton():
            crysSelect()
            parent.visual.togglePlane()

        def volSelect():
            parent.visual.setPlane(str(self.volDir.currentText()),self.volSel.value()-1.)

        def volButton():
            volSelect()
            parent.visual.togglePlane()

        def volDirection(i):
            self.volSel.setMaximum(self.shape[i])
            self.volSel.setTickInterval(self.shape[i]/10)
            self.volMax.setText(str(self.shape[i]))
            volSelect()

        vbox = QVBoxLayout()
        #crystal planes
        planeBut = QPushButton('Show/Hide')
        planeBut.clicked.connect(crysButton)
        hbox2 = QHBoxLayout()
        hbox2.addWidget(QLabel('h:'))
        hbox2.addWidget(QLabel('k:'))
        hbox2.addWidget(QLabel('l:'))
        hbox = QHBoxLayout()
        ph = QSpinBox()
        pk = QSpinBox()
        pl = QSpinBox()
        for i in [ph,pk,pl]:
            hbox.addWidget(i)
            i.valueChanged.connect(crysSelect)
            i.setMinimum(-99)
        vbox.addWidget(QLabel('Crystal plane:'))
        vbox.addLayout(hbox2)
        vbox.addLayout(hbox)
        vbox.addWidget(planeBut)
        vbox.addStretch()
        #volume slice
        volMin = QLabel('1')
        self.volDir = QComboBox()
        for i in 'xyz':
            self.volDir.addItem(i)
        self.volDir.currentIndexChanged.connect(volDirection)
        self.volDir.setDisabled(True)
        self.volSel = QSlider()
        self.volSel.setDisabled(True)
        self.volSel.setOrientation(1)
        self.volSel.setMinimum(1)
        self.volSel.setMaximum(1)
        self.volSel.setTickPosition(self.volSel.TicksBelow)
        self.volSel.setSingleStep(1)
        self.volSel.valueChanged.connect(volSelect)
        self.volMax = QLabel('1')
        self.volBut = QPushButton('Show/Hide')
        self.volBut.clicked.connect(volButton)
        self.volBut.setDisabled(True)
        self.shape = [1,1,1]
        hbox3=QHBoxLayout()
        hbox3.addWidget(volMin)
        hbox3.addWidget(self.volSel)
        hbox3.addWidget(self.volMax)
        hbox4=QHBoxLayout()
        hbox4.addWidget(QLabel('Plane:'))
        hbox4.addWidget(self.volDir)
        vbox.addWidget(QLabel('Volume slice:'))
        vbox.addLayout(hbox4)
        vbox.addLayout(hbox3)
        vbox.addWidget(self.volBut)
        vbox.addStretch()
        vbox.setContentsMargins(0,0,0,0)
        self.setLayout(vbox)

    def setMol(self,mol):
        vol = mol.get_vol()
        if vol is not None:
            self.shape = vol.shape
            lim=self.shape[self.volDir.currentIndex()]
            self.volSel.setMaximum(lim)
            self.volSel.setTickInterval(lim/10)
            self.volMax.setText(str(lim))
            self.volSel.setEnabled(True)
            self.volBut.setEnabled(True)
            self.volDir.setEnabled(True)
        else:
            self.volSel.setMaximum(1)
            self.volMax.setText('1')
            self.volSel.setDisabled(True)
            self.volBut.setDisabled(True)
            self.volDir.setDisabled(True)
