#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from collections import OrderedDict

from plane_tool import Plane
from vol_tool import Volume
from pick_tool import Picker
from script_tool import Script
from cellmod_tool import CellMod

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
                tools=OrderedDict([
                        ('Pick',Picker),
                        ('Script',Script),
                        ('Cell Mod.',CellMod),
                        ('Plane',Plane),
                        ('Volume',Volume)])
                for i in tools.items():
                    self.combo.addItem(i[0])
                    self.stack.addWidget(i[1](parent))

        def setMol(self,mol):
            if hasattr(self,'mol') and self.mol is mol:
                return
            self.mol = mol
            for i in range(self.stack.count()):
                self.stack.widget(i).setMol(mol)

        def setSel(self,sel):
            for i in range(self.stack.count()):
                self.stack.widget(i).setSel(sel)
