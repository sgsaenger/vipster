#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from numpy.linalg import norm
from numpy import degrees,arccos,dot,cross
from itertools import combinations

class ToolArea(QWidget):
        def __init__(self,parent):
                super(ToolArea,self).__init__()
                self.initStack()
                self.initPicker()
                self.initScript()
                self.initMult()
                #self.initPlane()
                self.parent = parent

        def setMol(self,mol):
                self.mol = mol
                #if hasattr(self.mol,'volume'):
                #        self.zPlane.setValidator(QIntValidator(0,len(self.mol.volume[0][0])))

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
            self.mol.set_bonds()
            self.parent.updateMolStep()


        def initPicker(self):
            self.pickArea = QTextEdit()
            self.pickArea.setReadOnly(True)
            tooltip = QLabel()
            tooltip.setText('Pick up to 4 atoms:')
            vbox=QVBoxLayout()
            vbox.addWidget(tooltip)
            vbox.addWidget(self.pickArea)
            pickWidget = QWidget()
            pickWidget.setLayout(vbox)
            self.combo.addItem('Pick')
            self.stack.addWidget(pickWidget)

        def pickHandler(self,sel):
            if len(sel)==0:
                self.pickArea.setPlainText('')
            elif len(sel)==1:
                at = self.mol.get_atom(sel[0],'angstrom')
                self.pickArea.setPlainText('Atom: '+str(sel[0])+'\n'+
                        'Type: '+at[0]+'\n'+
                        u'Coord(Å): {: 3.3f} {: 3.3f} {: 3.3f}'.format(*at[1][:]))
            else:
                at=[self.mol.get_atom(i,'angstrom') for i in sel]
                output='Atoms: '+str(sel)+'\nTypes: '
                for i in at:
                    output+=i[0]+' '
                diff12 = at[0][1]-at[1][1]
                output+=u'\nDist {1}-{2}:  {0:3.3f} Å'.format(norm(diff12),*sel)
                if len(sel)>2:
                    diff23 = at[1][1]-at[2][1]
                    output+=u'\nDist {1}-{2}:  {0:3.3f} Å'.format(norm(diff23),*sel[1:])
                if len(sel)>3:
                    diff34 = at[2][1]-at[3][1]
                    output+=u'\nDist {1}-{2}:  {0:3.3f} Å'.format(norm(diff34),*sel[2:])
                if len(sel)>2:
                    a123= degrees(arccos(dot(diff12,diff23)/(norm(diff12)*norm(diff23))))
                    output+=u'\nAngle {1}-{2}-{3}:  {0:3.3f}°'.format(a123,*sel)
                if len(sel)>3:
                    a234 = degrees(arccos(dot(diff23,diff34)/(norm(diff23)*norm(diff34))))
                    c123 = cross(diff12,diff23)
                    c234 = cross(diff23,diff34)
                    d1234 = degrees(arccos(dot(c123,c234)/(norm(c123)*norm(c234))))
                    output+=u'\nAngle {1}-{2}-{3}:  {0:3.3f}°'.format(a234,*sel[1:])
                    output+=u'\nDihedral {1}-{2}-{3}-{4}:  {0:3.3f}°'.format(d1234,*sel)
                self.pickArea.setPlainText(output)
            pass
