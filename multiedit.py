#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from numpy.linalg import norm
from numpy import degrees,arccos,dot,cross
from itertools import combinations

class ToolArea(QWidget):
        def __init__(self,parent):
                super(ToolArea,self).__init__()
                self.parent = parent
                self.stack = QStackedWidget()
                self.combo = QComboBox()
                self.combo.currentIndexChanged.connect(self.stack.setCurrentIndex)
                vbox = QVBoxLayout()
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('Tools:'))
                hbox.addWidget(self.combo)
                vbox.addLayout(hbox)
                vbox.addWidget(self.stack)
                self.setLayout(vbox)
                #initialize childwidgets (in order):
                self.initPicker()
                self.initScript()
                self.initMult()
                self.initVol()

        def setMol(self,mol):
                self.mol = mol
                #inform childwidgets about new molecule:
                self.volUpdate()
                self.pickUpdate()

####################################
# Multiply unit cell
####################################

        def initMult(self):
                self.mult = QWidget()
                self.combo.addItem('Mult. cell')
                self.stack.addWidget(self.mult)
                self.xmult = QLineEdit()
                self.xmult.setText('1')
                self.xmult.setValidator(QIntValidator(0,999))
                self.ymult = QLineEdit()
                self.ymult.setText('1')
                self.ymult.setValidator(QIntValidator(0,999))
                self.zmult = QLineEdit()
                self.zmult.setText('1')
                self.zmult.setValidator(QIntValidator(0,999))
                multBut = QPushButton('Apply')
                multBut.clicked.connect(self.multHandler)
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('x:'))
                hbox.addWidget(self.xmult)
                hbox.addWidget(QLabel('y:'))
                hbox.addWidget(self.ymult)
                hbox.addWidget(QLabel('z:'))
                hbox.addWidget(self.zmult)
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(multBut)
                self.mult.setLayout(vbox)

        def multHandler(self):
                if not hasattr(self,'mol'): return
                self.mol.mult(int(self.xmult.text()),int(self.ymult.text()),int(self.zmult.text()))
                self.mol.set_bonds()
                self.parent.updateMolStep()

####################################
# Volume-PP
####################################

        def initVol(self):
                self.vol = QWidget()
                self.combo.addItem('Volume')
                self.stack.addWidget(self.vol)
                volMin = QLabel('1')
                self.volSel = QSlider()
                self.volSel.setDisabled(True)
                self.volSel.setOrientation(1)
                self.volSel.setMinimum(1)
                self.volSel.setMaximum(1)
                self.volSel.setTickPosition(self.volSel.TicksBelow)
                self.volSel.setSingleStep(1)
                self.volMax = QLabel('1')
                self.volWarn = QLabel()
                self.volBut = QPushButton('Show/Hide')
                self.volBut.clicked.connect(self.volBHandler)
                self.volBut.setDisabled(True)
                hbox=QHBoxLayout()
                hbox.addWidget(volMin)
                hbox.addWidget(self.volSel)
                hbox.addWidget(self.volMax)
                vbox=QVBoxLayout()
                vbox.addWidget(QLabel('Plane:'))
                vbox.addLayout(hbox)
                vbox.addWidget(self.volBut)
                vbox.addWidget(self.volWarn)
                self.vol.setLayout(vbox)

        def volSHandler(self):
                self.parent.visual.setPlane('z',self.volSel.value()-1.)

        def volBHandler(self):
                self.volSHandler()
                self.parent.visual.togglePlane()

        def volUpdate(self):
                if hasattr(self.mol,'volume'):
                    lim=self.mol.get_vol().shape[2]
                    self.volSel.setMaximum(lim)
                    self.volSel.setTickInterval(lim/10)
                    self.volMax.setText(str(lim))
                    self.volSel.setEnabled(True)
                    self.volSel.valueChanged.connect(self.volSHandler)
                    self.volBut.setEnabled(True)
                else:
                    self.volSel.setMaximum(0)
                    self.volMax.setText('0')
                    self.volSel.valueChanged.disconnect()
                    self.volSel.setDisabled(True)
                    self.volBut.setDisabled(True)


####################################
# Script handling
####################################

        def initScript(self):
            scriptWidget = QWidget()
            self.combo.addItem('Script')
            self.stack.addWidget(scriptWidget)
            self.scriptArea = QTextEdit()
            self.scriptResult = QLabel()
            scriptBut = QPushButton('Evaluate')
            scriptBut.clicked.connect(self.scriptHandler)
            hbox=QHBoxLayout()
            hbox.addWidget(self.scriptResult)
            hbox.addWidget(scriptBut)
            vbox=QVBoxLayout()
            vbox.addWidget(self.scriptArea)
            vbox.addLayout(hbox)
            scriptWidget.setLayout(vbox)

        def scriptHandler(self):
            self.scriptResult.setText(self.mol.evalScript(str(self.scriptArea.toPlainText())))
            self.mol.set_bonds()
            self.parent.updateMolStep()

####################################
# Selected Atoms
####################################

        def initPicker(self):
            self.pickArea = QTextEdit()
            self.pickArea.setReadOnly(True)
            tooltip = QLabel()
            tooltip.setText('Pick up to 4 atoms:')
            self.pickWarn = QLabel()
            vbox=QVBoxLayout()
            vbox.addWidget(tooltip)
            vbox.addWidget(self.pickArea)
            vbox.addWidget(self.pickWarn)
            pickWidget = QWidget()
            pickWidget.setLayout(vbox)
            self.combo.addItem('Pick')
            self.stack.addWidget(pickWidget)

        def pickHandler(self,sel):
            br=0.52917721092
            self.pickWarn.setText('')
            if len(sel)==0:
                self.pickArea.setPlainText('')
            else:
                output ='Atoms: '+str([a[1] for a in sel])+'\n'
                output+='Types: '+str([a[2] for a in sel])+'\n'
                ids = [a[1] for a in sel]
                if len(sel)>1:
                    diff01 = sel[0][3]-sel[1][3]
                    output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff01)*br,*ids[:2])
                if len(sel)>2:
                    diff12 = sel[1][3]-sel[2][3]
                    output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff12)*br,*ids[1:3])
                    if len(sel)>3:
                        diff23 = sel[2][3]-sel[3][3]
                        output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff23)*br,*ids[2:])
                    a012 = degrees(arccos(dot(diff01,diff12)/(norm(diff01)*norm(diff12))))
                    output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a012,*ids[:3])
                if len(sel)>3:
                    a123 = degrees(arccos(dot(diff12,diff23)/(norm(diff12)*norm(diff23))))
                    c012 = cross(diff01,diff12)
                    c123 = cross(diff12,diff23)
                    d0123 = degrees(arccos(dot(c012,c123)/(norm(c012)*norm(c123))))
                    output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a123,*ids[1:])
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0:3.3f}°\n'.format(d0123,*ids)
                self.pickArea.setPlainText(output)

        def pickUpdate(self):
            if self.pickArea.toPlainText():
                self.pickWarn.setText('Data has changed!')
