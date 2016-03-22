# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt

from ..molecule import kpfmts
from .collapsiblewidget import collapsibleWidget

class MolTable(collapsibleWidget):
    def __init__(self,parent):
        self.parent = parent
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(['Type','x','y','z'])
        def showEvent(e):
            if hasattr(self,'mol'): self.fillTable
        self.table.showEvent = showEvent
        self.updatedisable = False
        def cellHandler(row,col):
            if self.updatedisable: return
            self.mol.initUndo()
            name = str(self.table.item(row,0).text())
            coord = [0,0,0]
            fix = [0,0,0]
            for j in [0,1,2]:
                coord[j]=float(self.table.item(row,j+1).text())
                fix[j]=int(self.table.item(row,j+1).checkState()/2)
            self.mol.setAtom(row,name,coord,fix=fix)
            self.mol.saveUndo('set coordinates')
            self.parent.updateMol()
        self.table.cellChanged.connect(cellHandler)
        def selHandler():
            if self.updatedisable: return
            sel = self.mol.getSelection()
            keep = []
            idx = set()
            for i in self.table.selectedRanges():
                for j in range(i.topRow(),i.bottomRow()+1):
                    idx.add(j)
            for i in idx:
                l = [k for k in sel if k[0] == i]
                if l:
                    keep.extend(l)
                else:
                    keep.append((i,(0,0,0)))
            self.mol.delSelection()
            for i in keep:
                self.mol.addSelection(i)
            self.parent.updateMol()
        self.table.itemSelectionChanged.connect(selHandler)
        super(MolTable,self).__init__('Coordinates',self.table)

    def setMol(self,mol):
        if self.sender() is self.table: return
        self.mol = mol
        if self.table.isVisible(): self.fillTable()

    def fillTable(self):
        self.updatedisable = True
        self.table.setRowCount(self.mol.nat)
        self.table.setVerticalHeaderLabels([str(x) for x in range(self.mol.nat)])
        for i in range(self.mol.nat):
            at = self.mol.getAtom(i)
            self.table.setItem(i,0,QTableWidgetItem(at[0]))
            for j in [0,1,2]:
                self.table.setItem(i,j+1,QTableWidgetItem(str(at[1][j])))
                self.table.item(i,j+1).setFlags(Qt.ItemFlag(51))
                self.table.item(i,j+1).setCheckState(int(at[3][j])*2)
        #update selection
        self.table.setSelectionMode(QAbstractItemView.MultiSelection)
        self.table.clearSelection()
        for j in set(i[0] for i in self.mol.getSelection()):
            self.table.selectRow(j)
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.updatedisable = False

class MolCell(collapsibleWidget):
    def __init__(self,parent):
        self.parent = parent
        self.container = QWidget()
        def showEvent(e):
            if hasattr(self,'mol'): self.fillContainer()
        self.container.showEvent = showEvent
        self.updatedisable = False
        super(MolCell,self).__init__('Cell Geometry',self.container,show=False)
        #Cell Dimension
        def cdmHandler():
            if self.updatedisable: return
            if float(self.cellDm.text()) == self.mol.getCellDim():return
            par = self.parent
            if self.appAll.isChecked():
                self.mol.setCellDimAll(float(self.cellDm.text()),self.scale.isChecked())
            else:
                self.mol.setCellDim(float(self.cellDm.text()),self.scale.isChecked())
            self.parent.updateMol()
        self.cellDm = QLineEdit()
        self.cellDm.setValidator(QDoubleValidator())
        self.cellDm.editingFinished.connect(cdmHandler)
        #Cell Vectors
        def vecHandler():
            if self.updatedisable: return
            vec=[[0,0,0],[0,0,0],[0,0,0]]
            for i in [0,1,2]:
                for j in [0,1,2]:
                    vec[i][j]=float(self.vtable.item(i,j).text())
            if vec == self.mol.getVec().tolist(): return
            par = self.parent
            if self.appAll.isChecked():
                self.mol.setVecAll(vec,self.scale.isChecked())
            else:
                self.mol.setVec(vec,self.scale.isChecked())
            self.parent.updateMol()
        self.vtable = QTableWidget()
        self.vtable.setColumnCount(3)
        self.vtable.setRowCount(3)
        self.vtable.setFixedHeight(120)
        self.vtable.setHorizontalHeaderLabels(['x','y','z'])
        self.vtable.itemChanged.connect(vecHandler)
        #Header Layout
        hbox = QHBoxLayout()
        hbox.addWidget(QLabel('Cell vectors:'))
        hbox.addStretch()
        hbox.addWidget(QLabel('Cell dimension:'))
        hbox.addWidget(self.cellDm)
        #Action modifiers
        self.appAll = QCheckBox('Apply to all Steps')
        self.scale = QCheckBox('Scale coordinates with cell')
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.appAll)
        hbox2.addWidget(self.scale)
        #Layout
        cellLayout = QVBoxLayout()
        cellLayout.addLayout(hbox)
        cellLayout.addWidget(self.vtable)
        cellLayout.addLayout(hbox2)
        self.container.setLayout(cellLayout)

    def setMol(self,mol):
        if self.sender() in [self.vtable,self.cellDm]: return
        self.mol = mol
        if self.container.isVisible(): self.fillContainer()

    def fillContainer(self):
        self.updatedisable = True
        self.cellDm.setText(str(self.mol.getCellDim()))
        vec = self.mol.getVec()
        for i in [0,1,2]:
            for j in [0,1,2]:
                self.vtable.setItem(i,j,QTableWidgetItem(str(vec[i,j])))
        self.updatedisable = False

class MolKPoints(collapsibleWidget):
    def __init__(self,parent):
        self.parent = parent
        self.container = QWidget()
        def showEvent(e):
            if hasattr(self,'mol'): self.fillContainer()
        self.container.showEvent = showEvent
        self.updatedisable = False
        super(MolKPoints,self).__init__('K-Points',self.container,show=False)
        #Monkhorst-Pack-Widget
        self.mpg = QWidget()
        self.mpg.widg = []
        def saveMPG():
            if self.updatedisable: return
            self.mol.setKpoints('mpg',[str(i.text()) for i in self.mpg.widg])
        for i in range(3):
            self.mpg.widg.append(QLineEdit())
            self.mpg.widg[-1].setValidator(QIntValidator(1,40000))
            self.mpg.widg[-1].setText('1')
            self.mpg.widg[-1].editingFinished.connect(saveMPG)
        for i in range(3):
            self.mpg.widg.append(QLineEdit())
            self.mpg.widg[-1].setValidator(QDoubleValidator())
            self.mpg.widg[-1].setText('0')
            self.mpg.widg[-1].editingFinished.connect(saveMPG)
        mpglayout = QGridLayout()
        mpglayout.addWidget(QLabel('x:'),0,0)
        mpglayout.addWidget(QLabel('y:'),0,2)
        mpglayout.addWidget(QLabel('z:'),0,4)
        mpglayout.addWidget(QLabel('x offs.:'),1,0)
        mpglayout.addWidget(QLabel('y offs.:'),1,2)
        mpglayout.addWidget(QLabel('z offs.:'),1,4)
        for i in range(6):
            mpglayout.addWidget(self.mpg.widg[i],i//3,2*(i%3)+1)
        self.mpg.setLayout(mpglayout)
        #Discrete-KPoints-Widget
        self.disc = QWidget()
        self.disc.table = QTableWidget()
        self.disc.table.setColumnCount(4)
        self.disc.table.setHorizontalHeaderLabels(['x','y','z','weight/count'])
        self.disc.table.setContextMenuPolicy(Qt.ActionsContextMenu)
        def cellHandler():
            if self.updatedisable: return
            disc = []
            for i in range(self.disc.table.rowCount()):
                disc.append([str(self.disc.table.item(i,j).text()) for j in range(4)])
            self.mol.setKpoints('discrete',disc)
        self.disc.table.cellChanged.connect(cellHandler)
        def newKpoint():
            self.updatedisable=True
            self.disc.table.setRowCount(self.disc.table.rowCount()+1)
            for i in range(4):
                self.disc.table.setItem(self.disc.table.rowCount()-1,i,QTableWidgetItem('1'))
            self.updatedisable=False
            cellHandler()
        newKp = QAction('New k-point',self.disc.table)
        newKp.setShortcut('Ctrl+K')
        newKp.triggered.connect(newKpoint)
        self.disc.table.addAction(newKp)
        def delKpoint():
            self.disc.table.removeRow(self.disc.table.currentRow())
            cellHandler()
        delKp = QAction('Delete k-point',self.disc.table)
        delKp.triggered.connect(delKpoint)
        self.disc.table.addAction(delKp)
        def optionHandler():
            self.mol.setKpoints('options',
                    {'bands':self.disc.bands.isChecked(),
                     'crystal':self.disc.crystal.isChecked()})
        self.disc.bands = QCheckBox()
        self.disc.bands.stateChanged.connect(optionHandler)
        self.disc.crystal = QCheckBox()
        self.disc.crystal.stateChanged.connect(optionHandler)
        disclayout = QGridLayout()
        disclayout.addWidget(self.disc.table,0,0,1,4)
        disclayout.addWidget(QLabel('Bands:'),1,0)
        disclayout.addWidget(self.disc.bands,1,1)
        disclayout.addWidget(QLabel('Crystal:'),1,2)
        disclayout.addWidget(self.disc.crystal,1,3)
        self.disc.setLayout(disclayout)
        #Format selector and widget-stack
        self.fmt = QComboBox()
        self.fmt.addItems(kpfmts)
        self.stack = QStackedWidget()
        self.stack.addWidget(QLabel('Gamma point only'))
        self.stack.addWidget(self.mpg)
        self.stack.addWidget(self.disc)
        self.fmt.currentIndexChanged.connect(self.changeKpoint)
        #Layout
        hbox = QHBoxLayout()
        hbox.addWidget(QLabel('Format:'))
        hbox.addWidget(self.fmt)
        vbox = QVBoxLayout()
        vbox.addLayout(hbox)
        vbox.addWidget(self.stack)
        self.container.setLayout(vbox)

    def changeKpoint(self,index):
        if self.updatedisable: return
        fmt = kpfmts[index]
        self.mol.setKpoints('active',fmt)
        self.stack.setCurrentIndex(index)
        if fmt == 'mpg':
            kp = self.mol.getKpoints(fmt)
            for i in range(6):
                self.mpg.widg[i].setText(kp[i])
        elif fmt == 'discrete':
            pass

    def setMol(self,mol):
        self.mol = mol
        if self.container.isVisible(): self.fillContainer()

    def fillContainer(self):
        self.updatedisable = True
        self.fmt.setCurrentIndex(kpfmts.index(self.mol.getKpoints('active')))
        self.updatedisable = False
        self.changeKpoint(self.fmt.currentIndex())
