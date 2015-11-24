# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from copy import copy
from ..ftypeplugins.lammpsData import lammps_atom_style,param as defaultParam

class LammpsDialog(QDialog):
    def __init__(self,parent=None,small=False):
        super(LammpsDialog,self).__init__(parent)
        self.lammps=LammpsTab(small)
        buttons=QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        vbox = QVBoxLayout()
        vbox.addWidget(self.lammps)
        vbox.addWidget(buttons)
        self.setLayout(vbox)

    @staticmethod
    def getWriteParams(parent=None,param=None):
        dialog = LammpsDialog(parent)
        if not param: param = defaultParam["default"]
        param = copy(param)
        dialog.lammps.setParam(param)
        dialog.exec_()
        return param

    @staticmethod
    def getReadParams(parent=None,param=None):
        dialog = LammpsDialog(parent,small)
        return param

class LammpsTab(QWidget):
    def __init__(self,small=False):
        super(LammpsTab,self).__init__()
        #self.fmt=QComboBox()
        #self.fmt.addItems(['angstrom','bohr'])
        #self.fmt.currentIndexChanged.connect(self.updateParam)
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
        grid.addWidget(QLabel("format"),0,0)
        grid.addWidget(self.fmt,0,1)
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
        self.param=param
        #self.atom.setCurrentIndex(["angstrom","bohr"].index(param["fmt"]))
        self.atom.setCurrentIndex(list(lammps_atom_style.keys()).index(param["atom_style"]))
        self.bond.setChecked(param["bonds"])
        self.angle.setChecked(param["angles"])
        self.dihedral.setChecked(param["dihedrals"])
        self.improper.setChecked(param["impropers"])

    def updateParam(self):
        #self.param["fmt"]=str(self.fmt.currentText())
        self.param["atom_style"]=str(self.atom.currentText())
        self.param["bonds"]=bool(self.bond.checkState())
        self.param["angles"]=bool(self.angle.checkState())
        self.param["dihedrals"]=bool(self.dihedral.checkState())
        self.param["impropers"]=bool(self.improper.checkState())
