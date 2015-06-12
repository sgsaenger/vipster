#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from tools import tools
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
                for i in tools.items():
                    self.combo.addItem(i[0])
                    self.stack.addWidget(i[1](parent))

        def setMol(self,mol):
            for i in range(self.stack.count()):
                self.stack.widget(i).setMol(mol)
