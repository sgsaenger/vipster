# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from ..ftypeplugins.lammpsData import lammps_atom_style

class LammpsTab(QWidget):
    def __init__(self):
        super(LammpsTab,self).__init__()
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
        grid.addWidget(QLabel("atom_style:"),0,0)
        grid.addWidget(self.atom,0,1)
        grid.addWidget(QLabel("Print bonds:"),1,0)
        grid.addWidget(self.bond,1,1)
        grid.addWidget(QLabel("Print angles:"),2,0)
        grid.addWidget(self.angle,2,1)
        grid.addWidget(QLabel("Print dihedrals:"),3,0)
        grid.addWidget(self.dihedral,3,1)
        grid.addWidget(QLabel("Print impropers:"),4,0)
        grid.addWidget(self.improper,4,1)
        self.setLayout(grid)

    def setParam(self,param):
        self.param=param
        self.atom.setCurrentIndex(list(lammps_atom_style.keys()).index(param["atom_style"]))
        self.bond.setChecked(param["bonds"])
        self.angle.setChecked(param["angles"])
        self.dihedral.setChecked(param["dihedrals"])
        self.improper.setChecked(param["impropers"])

    def updateParam(self):
        self.param["atom_style"]=str(self.atom.currentText())
        self.param["bonds"]=bool(self.bond.checkState())
        self.param["angles"]=bool(self.angle.checkState())
        self.param["dihedrals"]=bool(self.dihedral.checkState())
        self.param["impropers"]=bool(self.improper.checkState())
