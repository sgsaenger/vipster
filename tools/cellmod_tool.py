#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

####################################
# Modify unit cell
####################################

class CellMod(QWidget):

    def __init__(self,parent):
        super(CellMod,self).__init__()
        self.parent = parent
        vbox = QVBoxLayout()
        
        def modHandler():
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

        wrapBut = QPushButton('Wrap atoms')
        wrapBut.clicked.connect(modHandler)
        cropBut = QPushButton('Crop atoms')
        cropBut.clicked.connect(modHandler)
        vbox.addWidget(wrapBut)
        vbox.addWidget(cropBut)

        self.xmult = QSpinBox()
        self.xmult.setMinimum(1)
        self.ymult = QSpinBox()
        self.ymult.setMinimum(1)
        self.zmult = QSpinBox()
        self.zmult.setMinimum(1)
        mbox = QHBoxLayout()
        mbox.addWidget(QLabel('x:'))
        mbox.addWidget(self.xmult)
        mbox.addWidget(QLabel('y:'))
        mbox.addWidget(self.ymult)
        mbox.addWidget(QLabel('z:'))
        mbox.addWidget(self.zmult)
        vbox.addLayout(mbox)
        multBut = QPushButton('Multiply cell')
        multBut.clicked.connect(modHandler)
        vbox.addWidget(multBut)

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
        vbox.addWidget(self.reshape)
        rBut = QPushButton('Reshape cell')
        rBut.clicked.connect(modHandler)
        vbox.addWidget(rBut)

        vbox.addStretch()
        self.setLayout(vbox)

    def setMol(self,mol):
        self.mol = mol

    def setSel(self,sel):
        return
