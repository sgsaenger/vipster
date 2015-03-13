#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

class ToolArea(QWidget):
        def __init__(self,parent):
                super(ToolArea,self).__init__()
                self.initStack()
                self.initMult()
                self.initPlane()
                self.initScript()
                self.parent = parent

        def setMol(self,mol):
                self.mol = mol
                if hasattr(self.mol,'volume'):
                        self.zPlane.setValidator(QIntValidator(0,len(self.mol.volume[0][0])))

        def initStack(self):
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
                self.mol.set_pbc_bonds()
                self.parent.updateMolStep()

        def initPlane(self):
                self.plane = QWidget()
                self.combo.addItem('2D-PP')
                self.stack.addWidget(self.plane)
                self.zPlane = QLineEdit()
                self.zPlane.setText('0')
                self.zPlane.setValidator(QIntValidator(0,0))
                planeBut = QPushButton('Calc')
                planeBut.clicked.connect(self.planeHandler)
                hbox=QHBoxLayout()
                hbox.addWidget(QLabel('z-Plane:'))
                hbox.addWidget(self.zPlane)
                hbox.addWidget(planeBut)
                self.plane.setLayout(hbox)

        def planeHandler(self):
                #plane = self.mol.vol_plane(int(self.zPlane.text()))
                #print 'mean potential: '+str(self.mol.volume.mean())
                self.mol.stupid_me()
                #for k in range(self.mol.nvol[2]):
                        #plane = self.mol.vol_plane(k)
                        #img = QImage(len(plane),len(plane[0]),QImage.Format_RGB888)
                        #min=self.mol.volume.min()
                        #max=self.mol.volume.max()
                        #diff=max-min
                        #for i in range(len(plane)):
                        #        for j in range(len(plane[0])):
                        #                c=(plane[i][j]-min)/diff
                        #                r=255*np.exp(-6.25*(c-0.9)**2)
                        #                g=255*np.exp(-6.25*(c-0.5)**2)
                        #                b=255*np.exp(-6.25*(c-0.1)**2)
                        #                img.setPixel(i,j,qRgb(r,g,b))
                        #img.save('plane'+str(k)+'.png',None,-1)
                        #print "potential "+str(k)+': '+str(plane.mean())
                        #print "min: "+str(plane.min())+" max: "+str(plane.max())

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
