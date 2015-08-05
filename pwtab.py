#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt
from collections import OrderedDict

class PWTab(QTreeWidget):
    def __init__(self):
        super(PWTab,self).__init__()
        self.setColumnCount(2)
        self.setHeaderLabels(['Parameter','Value'])
        self.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.updatedisable=False
        self.itemChanged.connect(self.itemHandler)

        # Actions:
        self.addIons = QAction('Add Ions namelist',self)
        self.addIons.triggered.connect(self.createNamelist)
        self.addCell = QAction('Add Cell namelist',self)
        self.addCell.triggered.connect(self.createNamelist)
        self.addAction(self.addIons)
        self.addAction(self.addCell)
        newPar = QAction('New Parameter',self)
        newPar.setShortcut('Ctrl+N')
        newPar.triggered.connect(self.createParameter)
        self.addAction(newPar)
        delItem = QAction('Delete Item',self)
        delItem.setShortcut('Del')
        delItem.triggered.connect(self.deleteItem)
        self.addAction(delItem)

    def mouseDoubleClickEvent(self,e):
        #on double left click, edit selected item
        if (e.buttons() & 1):
            if self.currentItem().flags() & Qt.ItemIsEditable:
                self.editItem(self.currentItem(),self.currentColumn())

    def createNamelist(self):
        new = QTreeWidgetItem(self)
        if self.sender() is self.addIons:
            self.pw['&ions']=OrderedDict()
            new.setText(0,'&ions')
            self.addIons.setDisabled(True)
        elif self.sender() is self.addCell:
            self.pw['&cell']=OrderedDict()
            new.setText(0,'&cell')
            self.addCell.setDisabled(True)

    def createParameter(self):
        if self.currentItem().parent():
            new = QTreeWidgetItem(self.currentItem().parent())
        else:
            new = QTreeWidgetItem(self.currentItem())
        new.setFlags(new.flags()|Qt.ItemIsEditable)

    def deleteItem(self):
        if self.currentItem().parent():
            del self.pw[str(self.currentItem().parent().text(0))][str(self.currentItem().text(0))]
            self.currentItem().parent().removeChild(self.currentItem())
        elif self.currentItem().text(0) == '&ions':
            self.addIons.setEnabled(True)
            del self.pw['&ions']
            self.invisibleRootItem().removeChild(self.currentItem())
        elif self.currentItem().text(0) == '&cell':
            self.addCell.setEnabled(True)
            del self.pw['&cell']
            self.invisibleRootItem().removeChild(self.currentItem())

    def setPW(self,pw):
        self.pw=pw
        self.fillTree()

    def fillTree(self):
        self.updatedisable=True
        root = self.invisibleRootItem()
        #delete previous entries
        for i in range(root.childCount()):
            root.removeChild(root.child(0))
        #mandatory namelists
        for i in ['&control','&system','&electrons']:
            new = QTreeWidgetItem(self)
            new.setText(0,i)
        #show optional namelists only if existing
        if '&ions' in self.pw:
            new = QTreeWidgetItem(self)
            new.setText(0,'&ions')
            self.addIons.setDisabled(True)
        if '&cell' in self.pw:
            new = QTreeWidgetItem(self)
            new.setText(0,'&cell')
            self.addCell.setDisabled(True)
        #show child entries
        for i in range(root.childCount()):
            for j in self.pw[str(root.child(i).text(0))].items():
                new = QTreeWidgetItem(root.child(i))
                new.setText(0,j[0])
                new.setText(1,j[1])
                new.setFlags(new.flags()|Qt.ItemIsEditable)
        self.expandAll()
        self.resizeColumnToContents(0)
        self.updatedisable=False

    def itemHandler(self,item,column):
        if self.updatedisable:return
        if column:
            self.pw[str(item.parent().text(0))][str(item.text(0))]=str(item.text(1))
        elif item.parent():
            old=self.pw[str(item.parent().text(0))].keys()
            for i in range(item.parent().childCount()):
                if str(item.parent().child(i).text(0)) in old:
                    old.remove(str(item.parent().child(i).text(0)))
            self.pw[str(item.parent().text(0))]=\
                    OrderedDict([(str(item.text(0)),v) if k == old[0]\
                    else (k,v) for k,v in self.pw[str(item.parent().text(0))].items()])
