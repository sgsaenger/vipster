#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt

class PWTab(QSplitter):

        def __init__(self):
                super(PWTab,self).__init__()
                self.initTab()

        def initTab(self):
                self.tree = self.Tree()
                self.initKpoints()
                self.setOrientation(0)
                self.setChildrenCollapsible(False)
                self.setFixedWidth(445)
                self.addWidget(self.tree)
                self.addWidget(self.kp)
                self.setStretchFactor(0,1)

        #call fill functions when necessary:
        def setPW(self,pw):
                #save previous changes if applicable:
                if self.isVisible():
                        self.saveParam()
                #load parameters
                self.pw = pw
                #show if visible
                self.fillTree()
                self.fillKpoints()

        def showEvent(self,e):
                if hasattr(self,'pw'):
                        self.fillTree()
                        self.fillKpoints()

        def hideEvent(self,e):
                self.saveParam()

        def setCurrentIndex(self,i):
                if i>2: i=2
                self.kp.disp.setCurrentIndex(i)

        #init sections:
        def initKpoints(self):
                self.kp = QWidget()
                kp = self.kp
                #choose k point format
                kp.fmt = QComboBox()
                for i in ['gamma','automatic','tpiba','crystal','tpiba_b','crystal_b']:
                        kp.fmt.addItem(i)

                #Automatic: x,y,z, offset(x,y,z)
                kp.auto = QWidget()
                kp.auto.widg = []
                for i in [0,1,2]:
                        kp.auto.widg.append(QLineEdit())
                        kp.auto.widg[-1].setValidator(QIntValidator(1,40000))
                for i in [0,1,2]:
                        kp.auto.widg.append(QCheckBox())
                kp.auto.hbox1=QHBoxLayout()
                kp.auto.hbox1.addWidget(QLabel('x:'))
                kp.auto.hbox1.addWidget(kp.auto.widg[0])
                kp.auto.hbox1.addWidget(QLabel('y:'))
                kp.auto.hbox1.addWidget(kp.auto.widg[1])
                kp.auto.hbox1.addWidget(QLabel('z:'))
                kp.auto.hbox1.addWidget(kp.auto.widg[2])
                kp.auto.hbox2 = QHBoxLayout()
                kp.auto.hbox2.addWidget(QLabel('x offs.:'))
                kp.auto.hbox2.addWidget(kp.auto.widg[3])
                kp.auto.hbox2.addWidget(QLabel('y offs.:'))
                kp.auto.hbox2.addWidget(kp.auto.widg[4])
                kp.auto.hbox2.addWidget(QLabel('z offs.:'))
                kp.auto.hbox2.addWidget(kp.auto.widg[5])
                kp.auto.vbox = QVBoxLayout()
                kp.auto.vbox.addLayout(kp.auto.hbox1)
                kp.auto.vbox.addLayout(kp.auto.hbox2)
                kp.auto.setLayout(kp.auto.vbox)

                #discrete kpoints
                kp.disc = QTableWidget()
                kp.disc.setColumnCount(4)
                kp.disc.setHorizontalHeaderLabels(['x','y','z','weight'])
                kp.disc.setContextMenuPolicy(Qt.ActionsContextMenu)
                newKp = QAction('New k-point',kp.disc)
                newKp.setShortcut('Ctrl+K')
                newKp.triggered.connect(self.newKpoint)
                kp.disc.addAction(newKp)
                delKp = QAction('Delete k-point',kp.disc)
                delKp.triggered.connect(self.delKpoint)
                kp.disc.addAction(delKp)

                #stacked display of various formats
                kp.disp = QStackedWidget()
                kp.disp.addWidget(QLabel('Gamma point only'))
                kp.disp.addWidget(kp.auto)
                kp.disp.addWidget(kp.disc)
                kp.fmt.currentIndexChanged.connect(self.setCurrentIndex)

                #layout
                hbox = QHBoxLayout()
                hbox.addWidget(QLabel('K Points:'))
                hbox.addWidget(kp.fmt)
                vbox = QVBoxLayout()
                vbox.addLayout(hbox)
                vbox.addWidget(kp.disp)
                kp.setLayout(vbox)

        class Tree(QTreeWidget):
                def __init__(self):
                        super(PWTab.Tree,self).__init__()
                        self.setColumnCount(2)
                        self.setHeaderLabels(['Parameter','Value'])
                        self.setContextMenuPolicy(Qt.ActionsContextMenu)

                        # Actions:
                        newNl = QAction('New Namelist',self)
                        newNl.setShortcut('Ctrl+N')
                        newNl.triggered.connect(self.createNamelist)
                        self.addAction(newNl)
                        newPar = QAction('New Parameter',self)
                        newPar.setShortcut('Ctrl+P')
                        newPar.triggered.connect(self.createParameter)
                        self.addAction(newPar)
                        delItem = QAction('Delete Item',self)
                        delItem.setShortcut('Del')
                        delItem.triggered.connect(self.deleteItem)
                        self.addAction(delItem)

                def mouseDoubleClickEvent(self,e):
                        #on double left click, edit selected item
                        if (e.buttons() & 1):
                                self.editItem(self.currentItem(),self.currentColumn())

                def createNamelist(self):
                        new = QTreeWidgetItem(self)
                        new.setText(0,'&')
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                def createParameter(self):
                        #new = QTreeWidgetItem(self)
                        if self.currentItem().parent():
                                new = QTreeWidgetItem(self.currentItem().parent())
                        else:
                                new = QTreeWidgetItem(self.currentItem())
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                def deleteItem(self):
                        if self.currentItem().parent():
                                self.currentItem().parent().removeChild(self.currentItem())
                        else:
                                self.invisibleRootItem().removeChild(self.currentItem())

        #fill sections:
        def fillTree(self):
                root = self.tree.invisibleRootItem()
                #delete previous entries
                for i in range(root.childCount()):
                        root.removeChild(root.child(0))
                #mandatory namelists
                for i in ['&control','&system','&electrons']:
                        new = QTreeWidgetItem(self.tree)
                        new.setText(0,i)
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                #show optional namelists only if existing
                if '&ions' in self.pw:
                        new = QTreeWidgetItem(self.tree)
                        new.setText(0,'&ions')
                        new.setFlags(new.flags()|Qt.ItemIsEditable)
                if '&cell' in self.pw:
                        new = QTreeWidgetItem(self.tree)
                        new.setText(0,'&cell')
                        new.setFlags(new.flags()|Qt.ItemIsEditable)

                #show child entries
                for i in range(root.childCount()):
                        for j in self.pw[str(root.child(i).text(0))].items():
                                new = QTreeWidgetItem(root.child(i))
                                new.setText(0,j[0])
                                new.setText(1,j[1])
                                new.setFlags(new.flags()|Qt.ItemIsEditable)
                self.tree.expandAll()
                self.tree.resizeColumnToContents(0)

        def fillKpoints(self):
                self.kp.fmt.setCurrentIndex(['gamma','automatic','tpiba','crystal','tpiba_b','crystal_b'].index(self.pw['K_POINTS']['active']))
                if 'automatic' in self.pw['K_POINTS']:
                        for i in [0,1,2]:
                                self.kp.auto.widg[i].setText(self.pw['K_POINTS']['automatic'][i])
                        for i in [3,4,5]:
                                self.kp.auto.widg[i].setChecked(bool(int(self.pw['K_POINTS']['automatic'][i])))
                else:
                        for i in [0,1,2]:
                                self.kp.auto.widg[i].setText('')
                        for i in [3,4,5]:
                                self.kp.auto.widg[i].setChecked(False)
                if 'disc' in self.pw['K_POINTS']:
                        self.kp.disc.setRowCount(len(self.pw['K_POINTS']['disc']))
                        for i in range(len(self.pw['K_POINTS']['disc'])):
                                for j in [0,1,2,3]:
                                        self.kp.disc.setItem(i,j,QTableWidgetItem(self.pw['K_POINTS']['disc'][i][j]))
                else:
                        self.kp.disc.setRowCount(0)

        def newKpoint(self):
                self.kp.disc.setRowCount(self.kp.disc.rowCount()+1)

        def delKpoint(self):
                tab = self.kp.disc
                tab.removeRow(tab.currentRow())

        def saveParam(self):
                if not hasattr(self,'pw'):return
                #save monkhorst-pack-grid
                auto=[str(self.kp.auto.widg[i].text()) for i in [0,1,2]]
                auto+=[int(self.kp.auto.widg[i].isChecked()) for i in [3,4,5]]
                self.pw['K_POINTS']['automatic']=auto
                #save discrete k-points
                if self.kp.disc.rowCount() > 0:
                        disc = []
                        for i in range(self.kp.disc.rowCount()):
                                disc.append([self.kp.disc.item(i,j).text() for j in [0,1,2,3]])
                        self.pw['K_POINTS']['disc']=disc
                #save chosen k-point format
                self.pw['K_POINTS']['active']=str(self.kp.fmt.currentText())
                #save NameLists and parameters:
                for i in range(self.tree.invisibleRootItem().childCount()):
                        nl = self.tree.invisibleRootItem().child(i)
                        self.pw[str(nl.text(0))]={}
                        for i in range(nl.childCount()):
                                self.pw[str(nl.text(0))][str(nl.child(i).text(0))]=str(nl.child(i).text(1))
