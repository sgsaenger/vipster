# -*- coding: utf-8 -*-

from vipster.gui.qtwrapper import *

####################################
# Volume-PP
####################################


class Volume(QWidget):

    def __init__(self, parent):
        super(Volume, self).__init__()

        def setVal():
            if self.updatedisable:
                return
            if self.sender() is self.volVal:
                self.volSlide.blockSignals(True)
                self.volSlide.setValue(int(
                    (float(self.volVal.text()) - self.valRange.bottom()) /
                    (self.valRange.top() - self.valRange.bottom()) * 1000))
                self.volSlide.blockSignals(False)
            elif self.sender() is self.volSlide:
                self.volVal.setText(str(
                    self.volSlide.value() / 1000. *
                    (self.valRange.top() - self.valRange.bottom()) +
                    self.valRange.bottom()))
            elif self.sender() is self.volBut:
                self.showSurf = not self.showSurf
            parent.mol.setIsoVal(self.showSurf, float(self.volVal.text()),
                                 self.volBoth.isChecked())
            parent.updateMol()

        self.updatedisable = False
        self.showSurf = False
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
        self.volBut.setCheckable(True)
        self.volBut.clicked.connect(setVal)
        self.volVal = QLineEdit('0')
        self.volVal.setValidator(self.valRange)
        self.volVal.returnPressed.connect(setVal)
        self.volBoth = QCheckBox('+/-')
        self.volBoth.stateChanged.connect(setVal)
        hbox = QHBoxLayout()
        hbox.addWidget(self.volBoth)
        hbox.addWidget(self.volVal)
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.volMin)
        hbox2.addStretch()
        hbox2.addWidget(self.volMax)
        vbox = QVBoxLayout()
        vbox.addLayout(hbox)
        vbox.addWidget(self.volSlide)
        vbox.addLayout(hbox2)
        vbox.addWidget(self.volBut)
        vbox.addStretch()
        self.setLayout(vbox)

    def setMol(self, mol):
        vol = mol.getVol()
        self.updatedisable = True
        if vol is not None:
            self.valRange.setBottom(vol.min())
            self.valRange.setTop(vol.max())
            self.volMin.setText(str(vol.min()))
            self.volMax.setText(str(vol.max()))
            self.volSlide.setEnabled(True)
            self.volBut.setEnabled(True)
            self.volBoth.setEnabled(True)
            isoVal = mol.getIsoVal()
            self.volVal.setText(str(isoVal[1]))
            self.volSlide.setValue(int(
                (isoVal[1] - self.valRange.bottom()) /
                (self.valRange.top() - self.valRange.bottom()) * 1000))
            self.volBoth.setChecked(isoVal[2])
            self.volBut.setChecked(isoVal[0])
        else:
            self.volMin.setText('0')
            self.volMax.setText('0')
            self.volSlide.setDisabled(True)
            self.volSlide.setValue(0)
            self.volVal.setText('0')
            self.volBoth.setDisabled(True)
            self.volBut.setDisabled(True)
            self.volBoth.setChecked(False)
            self.volBut.setChecked(False)
        self.updatedisable = False
