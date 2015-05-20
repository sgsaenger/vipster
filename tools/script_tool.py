#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *

####################################
# Script handling
####################################

class Script(QWidget):
    def __init__(self,parent):
        super(Script,self).__init__()
        self.parent = parent

        def scriptHandler():
            self.scriptResult.setText(self.mol.evalScript(str(self.scriptArea.toPlainText())))
            self.mol.set_bonds()
            self.parent.updateMolStep()

        self.scriptArea = QTextEdit()
        self.scriptResult = QLabel()
        scriptBut = QPushButton('Evaluate')
        scriptBut.clicked.connect(scriptHandler)
        hbox=QHBoxLayout()
        hbox.addStretch()
        hbox.addWidget(scriptBut)
        vbox=QVBoxLayout()
        vbox.addWidget(self.scriptArea)
        vbox.addWidget(self.scriptResult)
        vbox.addLayout(hbox)
        self.setLayout(vbox)

    def setMol(self,mol):
        self.mol = mol

    def setSel(self,sel):
        return
