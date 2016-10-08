# -*- coding: utf-8 -*-

from vipster.gui.qtwrapper import *

####################################
# Script handling
####################################


class Script(QWidget):
    def __init__(self, parent):
        super(Script, self).__init__()
        self.parent = parent

        def scriptHandler():
            self.scriptResult.setText(
                self.mol.evalScript(str(self.scriptArea.toPlainText())))
            self.parent.updateMol()

        self.scriptArea = QTextEdit()
        self.scriptResult = QLabel()
        scriptBut = QPushButton('Evaluate')
        scriptBut.clicked.connect(scriptHandler)
        hbox = QHBoxLayout()
        hbox.addStretch()
        hbox.addWidget(scriptBut)
        vbox = QVBoxLayout()
        vbox.addWidget(self.scriptArea)
        vbox.addWidget(self.scriptResult)
        vbox.addLayout(hbox)
        vbox.setContentsMargins(0, 0, 0, 0)
        self.setLayout(vbox)
        self.setMinimumHeight(150)

    def setMol(self, mol):
        self.mol = mol
