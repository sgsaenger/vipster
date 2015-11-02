# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

####################################
# Volume-PP
####################################

class Volume(QWidget):

    def __init__(self,parent):
        super(Volume,self).__init__()

        def setVal():
            if self.sender() is self.volBoth:
                curVal=float(self.volVal.text())
            elif self.sender() is self.volVal:
                curVal=float(self.volVal.text())
                self.volSlide.blockSignals(True)
                self.volSlide.setValue(int((curVal-self.valRange.bottom())/(self.valRange.top()-self.valRange.bottom())*1000))
                self.volSlide.blockSignals(False)
            elif self.sender() is self.volSlide:
                curVal = self.volSlide.value()/1000.*(self.valRange.top()-self.valRange.bottom())+self.valRange.bottom()
                self.volVal.setText(str(curVal))
            parent.visual.setSurf(curVal,self.volBoth.isChecked())

        self.valRange = QDoubleValidator()
        self.valRange.setBottom(0.)
        self.valRange.setTop(0.)
        self.volMin = QLabel('0')
        self.volMax = QLabel('0')
        self.volSlide = QSlider()
        self.volSlide.setMinimum(0)
        self.volSlide.setMaximum(1000)
        self.volSlide.setSingleStep(1)
        self.volSlide.setOrientation(1)
        self.volSlide.setTickPosition(self.volSlide.TicksBelow)
        self.volSlide.setTickInterval(100)
        self.volSlide.valueChanged.connect(setVal)
        self.volBut = QPushButton('Show/Hide')
        self.volBut.clicked.connect(parent.visual.toggleSurf)
        self.volVal = QLineEdit('0')
        self.volVal.setValidator(self.valRange)
        self.volVal.returnPressed.connect(setVal)
        self.volBoth = QCheckBox('+/-')
        self.volBoth.stateChanged.connect(setVal)
        hbox=QHBoxLayout()
        hbox.addWidget(self.volBoth)
        hbox.addWidget(self.volVal)
        hbox2=QHBoxLayout()
        hbox2.addWidget(self.volMin)
        hbox2.addStretch()
        hbox2.addWidget(self.volMax)
        vbox=QVBoxLayout()
        vbox.addLayout(hbox)
        vbox.addWidget(self.volSlide)
        vbox.addLayout(hbox2)
        vbox.addWidget(self.volBut)
        self.setLayout(vbox)

    def setMol(self,mol):
        vol=mol.get_vol()
        if vol is not None:
            self.valRange.setBottom(vol.min())
            self.valRange.setTop(vol.max())
            self.volMin.setText(str(vol.min()))
            self.volMax.setText(str(vol.max()))
            self.volSlide.setEnabled(True)
            self.volBut.setEnabled(True)
        else:
            self.volMin.setText('0')
            self.volMax.setText('0')
            self.volSlide.setDisabled(True)
            self.volBut.setDisabled(True)
