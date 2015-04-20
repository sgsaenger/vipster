#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt
from copy import deepcopy

class MolArea(QWidget):

        def __init__(self,parent):
                super(MolArea,self).__init__()
                self.initTab()
                self.parent = parent

        def initTab(self):
                # coord fmt dropdown selector
                self.fmt = QComboBox()
                self.fmt.setToolTip('Select format of coordinates')
                for i in ['angstrom','bohr','crystal','alat']:
                        self.fmt.addItem(i)
                self.fmt.setCurrentIndex(2)
                self.fmt.currentIndexChanged.connect(self.fillTab)

                # layout1
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('Coordinates:'))
                hbox.addWidget(self.fmt)

                # show coordinates in table
                self.table = QTableWidget()
                self.table.setColumnCount(4)
                self.table.setHorizontalHeaderLabels(['Type','x','y','z'])
                self.table.cellChanged.connect(self.cellHandler)
                #show actions in right click menu
                self.table.setContextMenuPolicy(Qt.ActionsContextMenu)

                # show celldm
                self.cellDm = QLineEdit()
                self.cellDm.editingFinished.connect(self.cdmHandler)

                # show cell vectors in table
                self.vtable = QTableWidget()
                self.vtable.setColumnCount(3)
                self.vtable.setRowCount(3)
                self.vtable.setFixedHeight(120)
                self.vtable.setHorizontalHeaderLabels(['x','y','z'])
                self.vtable.itemChanged.connect(self.vecHandler)

                #cell header with label and celldm
                hbox2=QHBoxLayout()
                hbox2.addWidget(QLabel('Cell vectors:'))
                hbox2.addStretch(1)
                hbox2.addWidget(QLabel('Cell dimension:'))
                hbox2.addWidget(self.cellDm)

                #New Atom button:
                newB = QPushButton('New Atom')
                newB.clicked.connect(self.newAtom)

                #Copy Atom button:
                copyB = QPushButton('Copy Atom(s)')
                copyB.clicked.connect(self.copyAt)

                #paste button:
                pasteB = QPushButton('Paste Atom(s)')
                pasteB.clicked.connect(self.pasteAt)

                #delete button:
                delB = QPushButton('Delete Atom(s)')
                delB.clicked.connect(self.delAt)

                #sort buttons
                btns = QHBoxLayout()
                btns.addWidget(copyB)
                btns.addWidget(pasteB)
                btns.addWidget(newB)
                btns.addWidget(delB)

                # actions
                self.newA = QAction('New Atom',self)
                self.newA.setShortcut('Ctrl+N')
                self.newA.triggered.connect(self.newAtom)
                self.table.addAction(self.newA)
                self.copyA = QAction('Copy Atom(s)',self)
                self.copyA.setShortcut('Ctrl+C')
                self.copyA.triggered.connect(self.copyAt)
                self.table.addAction(self.copyA)
                self.pasteA = QAction('Paste Atom(s)',self)
                self.pasteA.setShortcut('Ctrl+V')
                self.pasteA.triggered.connect(self.pasteAt)
                self.table.addAction(self.pasteA)
                self.delA = QAction('Delete Atom(s)',self)
                self.delA.setShortcut('Del')
                self.delA.triggered.connect(self.delAt)
                self.table.addAction(self.delA)

                # Action modifiers
                self.appAll = QCheckBox('Apply to all Molecules')
                self.scale = QCheckBox('Scale coordinates with cell')
                hbox3 = QHBoxLayout()
                hbox3.addWidget(self.appAll)
                hbox3.addWidget(self.scale)

                # set Layout for Tab
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(self.table)
                vbox.addLayout(btns)
                vbox.addLayout(hbox2)
                vbox.addWidget(self.vtable)
                vbox.addLayout(hbox3)

                self.setLayout(vbox)
                self.resize(self.sizeHint())

        #load selected molecule
        def setMol(self,mol):
                #connect molecule
                self.mol = mol
                #fill if visible
                if self.isVisible(): self.fillTab()

        def showEvent(self,e):
                if hasattr(self,'mol'): self.fillTab()

        ##################################################################
        # EDIT FUNCTIONS
        ##################################################################
        def newAtom(self):
                self.mol.create_atom()
                self.mol.set_bonds()
                #update Main Widget
                self.parent.updateMolStep()

        def copyAt(self):
                self.sel = []
                for i in self.table.selectedRanges():
                        for j in range(i.topRow(),i.bottomRow()+1):
                                self.sel.append(self.mol.get_atom(j))

        def pasteAt(self):
                pos = self.table.currentRow()+1
                for at in reversed(self.sel):
                        self.mol.insert_atom(pos,at)
                self.mol.set_bonds()
                #update Main Widget
                self.parent.updateMolStep()

        def delAt(self):
                delrange = set()
                for i in self.table.selectedRanges():
                        for j in range(i.topRow(),i.bottomRow()+1):
                                delrange.add(j)
                for i in sorted(delrange,reverse=True):
                        self.mol.del_atom(i)
                self.mol.set_bonds()
                #update Main Widget
                self.parent.updateMolStep()

        #####################################################
        # UPDATE HANDLER
        ####################################################

        def cdmHandler(self):
                if self.updatedisable: return
                if float(self.cellDm.text()) == self.mol.get_celldm(self.fmt.currentText()):return
                par= self.parent
                if self.appAll.isChecked():
                    mols=[par.controller.get_mol(par.mlist.currentRow(),i) for i in range(int(par.maxStep.text()))]
                else:
                    mols=[self.mol]
                for m in mols:
                    m.set_celldm(float(self.cellDm.text()),self.scale.isChecked(),self.fmt.currentText())
                    m.set_bonds()

                #update Main Widget
                self.parent.updateMolStep()

        def cellHandler(self,row,col):
                if self.updatedisable: return
                #atom = self.table.currentRow()
                name = str(self.table.item(row,0).text())
                coord = [0,0,0]
                fix = [0,0,0]
                for j in [0,1,2]:
                        coord[j]=float(self.table.item(row,j+1).text())
                        fix[j]=int(not self.table.item(row,j+1).checkState()/2)
                        print(self.table.item(row,j+1).checkState()/2)
                print fix
                self.mol.set_atom(row,name,coord,self.fmt.currentText(),fix)
                print self.mol.get_atom(row)
                self.mol.set_bonds()

                #update Main Widget
                self.parent.updateMolStep()

        def vecHandler(self):
                if self.updatedisable: return
                vec=[[0,0,0],[0,0,0],[0,0,0]]
                for i in [0,1,2]:
                        for j in [0,1,2]:
                                vec[i][j]=float(self.vtable.item(i,j).text())
                if vec == self.mol.get_vec().tolist(): return
                par = self.parent
                if self.appAll.isChecked():
                    mols=[par.controller.get_mol(par.mlist.currentRow(),i) for i in range(int(par.maxStep.text()))]
                else:
                    mols=[self.mol]
                for m in mols:
                    m.set_vec(vec,self.scale.isChecked())
                    m.set_bonds()

                #update Main Widget
                self.parent.updateMolStep()

        ##############################################################
        # MAIN WIDGET UPDATE FUNCTION
        #############################################################

        def fillTab(self):
                #prevent handling of cell changes during fill
                self.updatedisable = True
                fmt = self.fmt.currentText()
                self.cellDm.setText(str(self.mol.get_celldm(fmt)))
                #fill atom table
                self.table.setRowCount(self.mol.get_nat())
                self.table.setVerticalHeaderLabels(map(str,range(self.mol.get_nat())))
                for i in range(self.mol.get_nat()):
                        at = self.mol.get_atom(i,fmt)
                        self.table.setItem(i,0,QTableWidgetItem(at[0]))
                        for j in [0,1,2]:
                                self.table.setItem(i,j+1,QTableWidgetItem(str(at[1][j])))
                                self.table.item(i,j+1).setFlags(Qt.ItemFlag(51))
                                self.table.item(i,j+1).setCheckState(int(not at[3][j])*2)
                #fill cell vec list
                vec = self.mol.get_vec()
                for i in [0,1,2]:
                        for j in [0,1,2]:
                                self.vtable.setItem(i,j,QTableWidgetItem(str(vec[i,j])))
                #reenable handling
                self.updatedisable = False

