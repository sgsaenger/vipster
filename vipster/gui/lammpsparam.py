# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from copy import copy
from ..ftypeplugins.lammpsData import lammps_atom_style,param as defaultParam

class LammpsParam(QWidget):
    def __init__(self,small=False):
        super(LammpsParam,self).__init__()
        self.updatedisable=False
        self.atom=QComboBox()
        self.atom.addItems(list(lammps_atom_style.keys()))
        self.atom.currentIndexChanged.connect(self.updateParam)
        self.bond=QCheckBox()
        self.bond.stateChanged.connect(self.updateParam)
        self.angle=QCheckBox()
        self.angle.stateChanged.connect(self.updateParam)
        self.dihedral=QCheckBox()
        self.dihedral.stateChanged.connect(self.updateParam)
        self.improper=QCheckBox()
        self.improper.stateChanged.connect(self.updateParam)
        grid = QGridLayout()
        grid.addWidget(QLabel("atom_style:"),1,0)
        grid.addWidget(self.atom,1,1)
        if not small:
            grid.addWidget(QLabel("Print bonds:"),2,0)
            grid.addWidget(self.bond,2,1)
            grid.addWidget(QLabel("Print angles:"),3,0)
            grid.addWidget(self.angle,3,1)
            grid.addWidget(QLabel("Print dihedrals:"),4,0)
            grid.addWidget(self.dihedral,4,1)
            grid.addWidget(QLabel("Print impropers:"),5,0)
            grid.addWidget(self.improper,5,1)
        self.setLayout(grid)

    def setParam(self,param):
        self.updatedisable=True
        self.param=param
        self.atom.setCurrentIndex(list(lammps_atom_style.keys()).index(param["atom_style"]))
        self.bond.setChecked(param["bonds"])
        self.angle.setChecked(param["angles"])
        self.dihedral.setChecked(param["dihedrals"])
        self.improper.setChecked(param["impropers"])
        self.updatedisable=False

    def updateParam(self):
        if self.updatedisable: return
        self.param["atom_style"]=str(self.atom.currentText())
        self.param["bonds"]=bool(self.bond.checkState())
        self.param["angles"]=bool(self.angle.checkState())
        self.param["dihedrals"]=bool(self.dihedral.checkState())
        self.param["impropers"]=bool(self.improper.checkState())
