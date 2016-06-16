# -*- coding: utf-8 -*-

from vipster.gui.qtwrapper import *

####################################
# Crystal and Volume planes
####################################

class Plane(QWidget):

    def __init__(self,parent):
        super(Plane,self).__init__()

        self.updatedisable = False
        vbox = QVBoxLayout()
        #crystal planes
        def setMilPlane():
            if self.updatedisable: return
            if self.sender() is self.milBut:
                self.showMilPlane = not self.showMilPlane
            parent.mol.setMilPlane(self.showMilPlane,[i.value() for i in self.milWidgets])
            parent.updateMol()
        self.showMilPlane = False
        vbox.addWidget(QLabel('Crystal plane:'))
        milGrid = QGridLayout()
        self.milWidgets = []
        for i,n in enumerate(['h','k','l']):
            self.milWidgets.append(QSpinBox())
            self.milWidgets[-1].valueChanged.connect(setMilPlane)
            milGrid.addWidget(QLabel(n),0,i)
            milGrid.addWidget(self.milWidgets[-1],1,i)
        vbox.addLayout(milGrid)
        self.milBut = QPushButton('Show/Hide')
        self.milBut.clicked.connect(setMilPlane)
        self.milBut.setCheckable(True)
        vbox.addWidget(self.milBut)
        vbox.addStretch()
        #volume slice
        def setVolPlane():
            if self.updatedisable: return
            if self.sender() is self.volDir:
                self.volSel.blockSignals(True)
                d = self.volDir.currentIndex()
                self.volSel.setMaximum(self.shape[d])
                self.volSel.setTickInterval(self.shape[d]/10)
                self.volMax.setText(str(self.shape[d]))
                self.volSel.blockSignals(False)
            elif self.sender() is self.volBut:
                self.showVolPlane = not self.showVolPlane
            parent.mol.setVolPlane(self.showVolPlane,self.volDir.currentIndex(),float(self.volSel.value()-1))
            parent.updateMol()
        self.showVolPlane = False
        volMin = QLabel('1')
        self.volDir = QComboBox()
        for i in 'xyz':
            self.volDir.addItem(i)
        self.volDir.currentIndexChanged.connect(setVolPlane)
        self.volDir.setDisabled(True)
        self.volSel = QSlider()
        self.volSel.setDisabled(True)
        self.volSel.setOrientation(1)
        self.volSel.setMinimum(1)
        self.volSel.setMaximum(1)
        self.volSel.setTickPosition(self.volSel.TicksBelow)
        self.volSel.setSingleStep(1)
        self.volSel.valueChanged.connect(setVolPlane)
        self.volMax = QLabel('1')
        self.volBut = QPushButton('Show/Hide')
        self.volBut.clicked.connect(setVolPlane)
        self.volBut.setCheckable(True)
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
        self.updatedisable = True
        vol = mol.getVol()
        milPlane = mol.getMilPlane()
        self.milBut.setChecked(milPlane[0])
        for i in range(3):
            self.milWidgets[i].setValue(milPlane[1][i])
        if vol is not None:
            self.shape = vol.shape
            lim=self.shape[self.volDir.currentIndex()]
            volPlane = mol.getVolPlane()
            self.volSel.setMaximum(lim)
            self.volSel.setTickInterval(lim/10)
            self.volSel.setEnabled(True)
            self.volSel.setValue(volPlane[2]+1)
            self.volMax.setText(str(lim))
            self.volBut.setEnabled(True)
            self.volBut.setChecked(volPlane[0])
            self.volDir.setEnabled(True)
            self.volDir.setCurrentIndex(volPlane[1])
        else:
            self.volSel.setValue(1)
            self.volSel.setDisabled(True)
            self.volMax.setText('1')
            self.volDir.setCurrentIndex(0)
            self.volDir.setDisabled(True)
            self.volBut.setDisabled(True)
            self.volBut.setChecked(False)
        self.updatedisable = False
