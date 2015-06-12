#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from PyQt4.QtCore import Qt

class PWTab(QTreeWidget):
    def __init__(self):
        super(PWTab,self).__init__()
        self.setColumnCount(2)
        self.setHeaderLabels(['Parameter','Value'])
        self.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.updatedisable=False
        self.itemChanged.connect(self.itemHandler)

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

    def setPW(self,pw):
        self.updatedisable=True
        self.pw=pw
        self.fillTree()
        self.updatedisable=False

    def fillTree(self):
        root = self.invisibleRootItem()
        #delete previous entries
        for i in range(root.childCount()):
            root.removeChild(root.child(0))
        #mandatory namelists
        for i in ['&control','&system','&electrons']:
            new = QTreeWidgetItem(self)
            new.setText(0,i)
            new.setFlags(new.flags()|Qt.ItemIsEditable)

        #show optional namelists only if existing
        if '&ions' in self.pw:
            new = QTreeWidgetItem(self)
            new.setText(0,'&ions')
            new.setFlags(new.flags()|Qt.ItemIsEditable)
        if '&cell' in self.pw:
            new = QTreeWidgetItem(self)
            new.setText(0,'&cell')
            new.setFlags(new.flags()|Qt.ItemIsEditable)

        #show child entries
        for i in range(root.childCount()):
            for j in self.pw[str(root.child(i).text(0))].items():
                new = QTreeWidgetItem(root.child(i))
                new.setText(0,j[0])
                new.setText(1,j[1])
                new.setFlags(new.flags()|Qt.ItemIsEditable)
        self.expandAll()
        self.resizeColumnToContents(0)

    def itemHandler(self,item):
        if self.updatedisable:return
        if item.parent() == None:
            self.pw[str(item.text(0))]={}
            for i in range(item.childCount()):
                self.pw[str(item.text(0))][str(item.child(i).text(0))]=str(item.child(i).text(1))
        else:
            self.pw[str(item.parent().text(0))][str(item.text(0))]=str(item.text(1))
