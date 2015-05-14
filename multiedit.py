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
                #self.initPicker()
                #self.initScript()
                #self.initMod()
                tools={'Plane':self.Plane,
                        'Volume':self.Volume}
                for i in tools.items():
                    self.combo.addItem(i[0])
                    self.stack.addWidget(i[1](parent))

        def setMol(self,mol):
            for i in range(self.stack.count()):
                self.stack.widget(i)._update(mol)
                #self.mol = mol
                #inform childwidgets about new molecule:
                #self.volUpdate()
                #self.planeUpdate()
                #self.pickUpdate()

####################################
# Modify unit cell
####################################

        def initMod(self):
            self.mult = QWidget()
            self.combo.addItem('Mod. cell')
            self.stack.addWidget(self.mult)

            wrapBut = QPushButton('Wrap atoms')
            wrapBut.clicked.connect(self.modHandler)
            cropBut = QPushButton('Crop atoms')
            cropBut.clicked.connect(self.modHandler)

            self.xmult = QSpinBox()
            self.xmult.setMinimum(1)
            self.ymult = QSpinBox()
            self.ymult.setMinimum(1)
            self.zmult = QSpinBox()
            self.zmult.setMinimum(1)
            multBut = QPushButton('Multiply cell')
            multBut.clicked.connect(self.modHandler)
            mbox = QHBoxLayout()
            mbox.addWidget(QLabel('x:'))
            mbox.addWidget(self.xmult)
            mbox.addWidget(QLabel('y:'))
            mbox.addWidget(self.ymult)
            mbox.addWidget(QLabel('z:'))
            mbox.addWidget(self.zmult)

            self.reshape = QTableWidget()
            self.reshape.setColumnCount(3)
            self.reshape.setRowCount(3)
            self.reshape.setFixedHeight(120)
            self.reshape.setColumnWidth(0,84)
            self.reshape.setColumnWidth(1,84)
            self.reshape.setColumnWidth(2,84)
            self.reshape.setHorizontalHeaderLabels(['x','y','z'])
            for i in range(3):
                for j in range(3):
                    self.reshape.setItem(i,j,QTableWidgetItem(str(0.0)))
            rBut = QPushButton('Reshape cell')
            rBut.clicked.connect(self.modHandler)

            vbox = QVBoxLayout()
            vbox.addWidget(wrapBut)
            vbox.addWidget(cropBut)
            vbox.addLayout(mbox)
            vbox.addWidget(multBut)
            vbox.addWidget(self.reshape)
            vbox.addWidget(rBut)
            vbox.addStretch()
            self.mult.setLayout(vbox)

        def modHandler(self):
            if not hasattr(self,'mol'): return
            reason=self.sender().text()
            if reason=='Multiply cell':
                self.mol.mult(int(self.xmult.text()),int(self.ymult.text()),int(self.zmult.text()))
            elif reason=='Reshape cell':
                vec=[[0,0,0],[0,0,0],[0,0,0]]
                for i in [0,1,2]:
                        for j in [0,1,2]:
                                vec[i][j]=float(self.reshape.item(i,j).text())
                self.mol.reshape(vec)
            elif reason=='Wrap atoms':
                self.mol.wrap()
            elif reason=='Crop atoms':
                self.mol.crop()
            self.mol.set_bonds()
            self.parent.updateMolStep()

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
                    diff12 = sel[2][3]-sel[1][3]
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

####################################
# Crystal and Volume planes
####################################

        class Plane(QWidget):

            def __init__(self,parent):
                super(ToolArea.Plane,self).__init__()

                def crysSHandler():
                    parent.visual.setPlane('c',[
                        1./ph.value() if ph.value() else 0,
                        1./pk.value() if pk.value() else 0,
                        1./pl.value() if pl.value() else 0])

                def crysBHandler():
                    crysSHandler()
                    parent.visual.togglePlane()

                def volSHandler():
                    parent.visual.setPlane(str(self.volDir.currentText()),self.volSel.value()-1.)

                def volBHandler():
                    volSHandler()
                    parent.visual.togglePlane()

                vbox = QVBoxLayout()
                #crystal planes
                planeBut = QPushButton('Show/Hide')
                planeBut.clicked.connect(crysBHandler)
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
                    i.valueChanged.connect(crysSHandler)
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
                self.volDir.currentIndexChanged.connect(volSHandler)
                self.volDir.setDisabled(True)
                self.volSel = QSlider()
                self.volSel.setDisabled(True)
                self.volSel.setOrientation(1)
                self.volSel.setMinimum(1)
                self.volSel.setMaximum(1)
                self.volSel.setTickPosition(self.volSel.TicksBelow)
                self.volSel.setSingleStep(1)
                self.volSel.valueChanged.connect(volSHandler)
                self.volMax = QLabel('1')
                self.volBut = QPushButton('Show/Hide')
                self.volBut.clicked.connect(volBHandler)
                self.volBut.setDisabled(True)
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
                self.setLayout(vbox)

            def _update(self,mol):
                if hasattr(mol,'_vol'):
                    lim=mol.get_vol().shape[ord(str(self.volDir.currentText()))-120]
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

####################################
# Volume-PP
####################################

        class Volume(QWidget):

            def __init__(self,parent):
                super(ToolArea.Volume,self).__init__()

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
                    pass

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

            def _update(self,mol):
                if hasattr(mol,'_vol'):
                    self.vol = mol.get_vol()
                    self.volMin.setText(str(self.vol.min()))
                    self.volSel.setMaximum(1000)
                    self.volMax.setText(str(self.vol.max()))
                    self.volSel.setEnabled(True)
                    self.volBut.setEnabled(True)
                    self.volSel.setValue(0)
                else:
                    self.volMin.setText('0')
                    self.volSel.setMaximum(0)
                    self.volMax.setText('0')
                    self.volSel.setDisabled(True)
                    self.volBut.setDisabled(True)
