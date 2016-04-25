# -*- coding: utf-8 -*-

from os.path import splitext
from os import getcwd
from copy import deepcopy

from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer,Qt,QSettings,QByteArray
try:
    from PyQt4.QtCore import QString
except:
    QString = str

from vipster.gui.viewport import VisualWidget
from vipster.gui.tools import tools
from vipster.gui.parameterwidgets import ParamTab,ParamDialog
from vipster.gui.molwidgets import MolCell,MolKPoints,MolTable
from vipster.gui.confwidgets import Settings,PseGlobal,PseMol

from vipster.molecule import fmts,Molecule
from vipster.settings import saveConfig
from vipster.iowrapper import availParam,newParam,_guiInNames,_guiOutNames,_paramdict,readFile,writeFile

GuiMolecules = []
GuiParameters = []

def launchVipster(m=GuiMolecules,p=GuiParameters):
    app = QApplication([])
    gui = MainWindow(m,p)
    app.aboutToQuit.connect(gui.deleteLater)
    app.exec_()

__all__=['GuiMolecules','GuiParameters','launchVipster']

class MainWindow(QMainWindow):

    def __init__(self,m,p):
        super(MainWindow,self).__init__()
        self.setWindowTitle("Vipster")
        #
        self.molecules = m
        if not self.molecules:
            self.molecules.append(Molecule())
        self.parameters = p
        self.mol = self.molecules[0]
        self.param = None
        #
        self.visual = VisualWidget(self)
        self.setCentralWidget(self.visual)
        #
        self.initEditWidgets()
        #
        self.initTools()
        #
        self.initFmt()
        #
        self.initFileMenu()
        #
        self.initEditMenu()
        #
        self.initFinish()

    def initEditWidgets(self):
        #Molecule
        self.mlist = QListWidget()
        self.mlist.currentRowChanged.connect(self.updateMol)
        self.mtable = MolTable(self)
        self.mcell = MolCell(self)
        self.mkpoints = MolKPoints(self)
        mscroll = QScrollArea()
        mwidget = QWidget()
        mlayout = QVBoxLayout()
        mlayout.addWidget(self.mlist)
        mlayout.addWidget(self.mtable)
        mlayout.addWidget(self.mcell)
        mlayout.addWidget(self.mkpoints)
        mwidget.setLayout(mlayout)
        mscroll.setWidget(mwidget)
        mscroll.setWidgetResizable(True)
        mscroll.setFrameStyle(0)
        mdock = QDockWidget("Molecules")
        mdock.setWidget(mscroll)
        mdock.setObjectName("Molecules")
        self.addDockWidget(1,mdock)
        #Parameter
        self.plist = QListWidget()
        self.plist.currentRowChanged.connect(self.updateParam)
        self.pedit = ParamTab()
        pwidget = QWidget()
        playout = QVBoxLayout()
        playout.addWidget(self.plist)
        playout.addWidget(self.pedit)
        pwidget.setLayout(playout)
        pdock = QDockWidget("Parameters")
        pdock.setWidget(pwidget)
        pdock.setObjectName("Parameters")
        self.addDockWidget(1,pdock)
        #Settings and PSE
        cdock = QDockWidget("Settings")
        cdock.setWidget(Settings(self))
        cdock.setObjectName("Settings")
        self.addDockWidget(1,cdock)
        self.pseglob = PseGlobal(self)
        pgdock = QDockWidget("PSE")
        pgdock.setWidget(self.pseglob)
        pgdock.setObjectName("PSE")
        self.addDockWidget(1,pgdock)
        self.psemol = PseMol(self)
        pmdock = QDockWidget("PSE (Mol.)")
        pmdock.setWidget(self.psemol)
        pmdock.setObjectName("PSE (Mol.)")
        self.addDockWidget(1,pmdock)
        self.setTabPosition(Qt.AllDockWidgetAreas,QTabWidget.North)
        self.tabifyDockWidget(cdock,pgdock)
        self.tabifyDockWidget(pgdock,pmdock)
        self.tabifyDockWidget(pmdock,pdock)
        self.tabifyDockWidget(pdock,mdock)

    def initTools(self):
        self.tools = []
        tb = self.addToolBar("Tools")
        tb.setObjectName("Tools")
        tb.addWidget(QLabel("Tools:"))
        for i in tools.items():
            dw = QDockWidget(i[0])
            dw.setWidget(i[1](self))
            dw.setObjectName(i[0])
            self.addDockWidget(2,dw)
            dw.hide()
            tb.addAction(dw.toggleViewAction())
            self.tools.append(dw.widget())

    def initFmt(self):
        self.fmt = QComboBox()
        self.fmt.addItems(fmts)
        self.fmt.currentIndexChanged[QString].connect(self.updateFmt)
        fb = self.addToolBar("Format")
        fb.setObjectName("Format")
        fb.addWidget(self.fmt)

    def initFileMenu(self):
        fMenu = self.menuBar().addMenu("&File")
        #
        newAction = QAction("&New Molecule",self)
        newAction.setShortcut("Ctrl+N")
        newAction.triggered.connect(self.newMol)
        fMenu.addAction(newAction)
        #
        nMenu = fMenu.addMenu("From existing Molecule:")
        copyStep = QAction("&Copy single Step",self)
        copyStep.triggered.connect(self.newMol)
        nMenu.addAction(copyStep)
        copyTraj = QAction("&Copy Trajectory",self)
        copyTraj.triggered.connect(self.newMol)
        nMenu.addAction(copyTraj)
        #
        loadAction = QAction("&Load Molecule(s)",self)
        loadAction.setShortcut("Ctrl+O")
        loadAction.triggered.connect(self.loadMol)
        fMenu.addAction(loadAction)
        #
        saveAction = QAction("&Save Molecule",self)
        saveAction.setShortcut("Ctrl+S")
        saveAction.triggered.connect(self.saveMol)
        fMenu.addAction(saveAction)
        #
        fMenu.addSeparator()
        pMenu = fMenu.addMenu("New &Parameter set")
        #check for available parameter sets and add
        for guiname,cliname in _guiOutNames.items():
            p = availParam(cliname)
            if not p:
                continue
            if len(p)==1:
                action=QAction(guiname,self)
                action.triggered.connect(self.newParam)
                pMenu.addAction(action)
            else:
                p2Menu = pMenu.addMenu(guiname)
                for i in p:
                    action=QAction(i,p2Menu)
                    action.triggered.connect(self.newParam)
                    p2Menu.addAction(action)
        #
        saveParamAction = QAction("Save Parameter set",self)
        saveParamAction.triggered.connect(self.saveParam)
        fMenu.addAction(saveParamAction)
        #Screenshot
        fMenu.addSeparator()
        scrotAction = QAction("&Screenshot",self)
        scrotAction.setShortcut("Ctrl+P")
        scrotAction.triggered.connect(self.savePic)
        fMenu.addAction(scrotAction)
        #
        fMenu.addSeparator()
        exitAction = QAction("&Exit",self)
        exitAction.setShortcut("Ctrl+Q")
        exitAction.triggered.connect(self.close)
        fMenu.addAction(exitAction)

    def initEditMenu(self):
        eMenu = self.menuBar().addMenu("&Edit")
        #
        def newFun():
            self.mol.initUndo()
            self.mol.newAtom()
            self.mol.saveUndo('create atom')
            self.updateMol()
        newAction = QAction('&New atom',self)
        newAction.setShortcut('n')
        newAction.triggered.connect(newFun)
        eMenu.addAction(newAction)
        #
        def delFun():
            sel = set(i[0] for i in self.mol.getSelection())
            if not sel: return
            self.mol.initUndo()
            for i in sorted(sel,reverse=True):
                self.mol.delAtom(i)
            self.mol.saveUndo('delete atom(s)')
            self.updateMol()
        delAction = QAction('&Delete atom(s)',self)
        delAction.setShortcut(QKeySequence.Delete)
        delAction.triggered.connect(delFun)
        eMenu.addAction(delAction)
        #
        def copyFun():
            sel = set(i[0] for i in self.mol.getSelection())
            self.copyBuf = []
            for i in sel:
                self.copyBuf.append(deepcopy(self.mol.getAtom(i,fmt='bohr')))
        copyAction = QAction('&Copy atom(s)',self)
        copyAction.setShortcut('Ctrl+C')
        copyAction.triggered.connect(copyFun)
        eMenu.addAction(copyAction)
        #
        def pasteFun():
            if not self.copyBuf: return
            self.mol.initUndo()
            for i in self.copyBuf:
                self.mol.newAtom(*i,fmt='bohr')
            self.mol.saveUndo('copy atom(s)')
            self.updateMol()
        pasteAction = QAction('&Paste atom(s)',self)
        pasteAction.setShortcut('Ctrl+V')
        pasteAction.triggered.connect(pasteFun)
        eMenu.addAction(pasteAction)
        #
        eMenu.addSeparator()
        def hideFun():
            sel = set(i[0] for i in self.mol.getSelection())
            for i in sel:
                self.mol.setAtom(i,hidden=True)
            self.updateMol()
        hideAction = QAction('&Hide atom(s)',self)
        hideAction.setShortcut('h')
        hideAction.triggered.connect(hideFun)
        eMenu.addAction(hideAction)
        def showFun():
            sel = set(i[0] for i in self.mol.getSelection())
            for i in sel:
                self.mol.setAtom(i,hidden=False)
            self.updateMol()
        showAction = QAction('&Show atom(s)',self)
        showAction.setShortcut('j')
        showAction.triggered.connect(showFun)
        eMenu.addAction(showAction)
        #
        eMenu.addSeparator()
        def undoFun():
            self.mol.undo()
            self.updateMol()
        self.undoAction = QAction('Undo',self)
        self.undoAction.setShortcut('Ctrl+Z')
        self.undoAction.triggered.connect(undoFun)
        eMenu.addAction(self.undoAction)

    def initFinish(self):
        for m in self.molecules:
            self.mlist.addItem(m.name)
        for p in self.parameters:
            self.plist.addItem(p["name"])
        self.show()
        self.mlist.setCurrentRow(self.mlist.count()-1)
        self.plist.setCurrentRow(self.plist.count()-1)
        settings = QSettings("hein09","vipster")
        g = settings.value("geometry")
        w = settings.value("windowState")
        if g:
            self.restoreGeometry(g if isinstance(g,QByteArray) else g.toByteArray())
        if g:
            self.restoreState(w if isinstance(w,QByteArray) else w.toByteArray())

### I/O FUNCTIONS ###

    def newMol(self):
        sender = self.sender().text()
        if "Molecule" in sender:
            mol = Molecule()
        elif "Step" in sender:
            mol = self.mol.copy()
        elif "Trajectory" in sender:
            mol = self.mol.copyAll()
        self.molecules.append(mol)
        self.mlist.addItem(mol.name)
        self.mlist.setCurrentRow(self.mlist.count()-1)

    def loadMol(self):
        fname = QFileDialog.getOpenFileName(self,"Open File",getcwd())
        if not fname: return
        ftype = QInputDialog.getItem(self,"Choose File type","File type:",list(_guiInNames.keys()),0,False)
        if not ftype[1]: return
        ftype = _guiInNames[str(ftype[0])]
        m,p = readFile(fname,ftype)
        self.molecules.append(m)
        self.mlist.addItem(m.name)
        self.mlist.setCurrentRow(self.mlist.count()-1)
        if p:
            self.parameters.append(p)
            self.plist.addItem(p["name"])
            self.plist.setCurrentRow(self.plist.count()-1)

    def saveMol(self):
        fname = QFileDialog.getSaveFileName(self,"Save File",getcwd())
        if not fname: return
        ftype = QInputDialog.getItem(self,"Choose File type","File type:",list(_guiOutNames.keys()),0,False)
        if not ftype[1]: return
        ftype = _guiOutNames[str(ftype[0])]
        try:
            writeFile(self.mol,ftype,fname,self.param)
        except Exception as error:
            pd = ParamDialog(ftype,self.parameters)
            if pd.exec_():
                writeFile(self.mol,ftype,fname,pd.getParam())

    def newParam(self):
        if self.sender().parent() != self:
            prog = _guiOutNames[str(self.sender().parent().title())]
            var = str(self.sender().text())
        else:
            prog = _guiOutNames[str(self.sender().text())]
            var = "default"
        param = newParam(prog,var)
        self.parameters.append(param)
        self.plist.addItem(param["name"])
        self.plist.setCurrentRow(self.plist.count()-1)

    def saveParam(self):
        if self.param:
            p = self.param
            t = p["type"]
            n = QInputDialog.getText(self,"Choose Name","Name:")
            if not n[1]: return
            n = str(n[0])
            _paramdict[t][n]=p
            _paramdict[t][n]["name"]=n
            saveConfig()

    def savePic(self):
        fn = QFileDialog.getSaveFileName(self,"Save Screenshot",getcwd(),"Portable Network Graphics (*.png)")
        if not fn: return
        if splitext(str(fn))[1] == "": fn+=".png"
        img = self.visual.makePicture()
        img.save(fn,None,0)

    def closeEvent(self,e):
        settings = QSettings("hein09","vipster")
        settings.setValue("geometry",self.saveGeometry())
        settings.setValue("windowState",self.saveState())
        super(QMainWindow,self).closeEvent(e)

### PUBLIC UPDATE FUNCTIONS ###

    def updateFmt(self,fmt):
        self.mol.setFmt(fmt)
        self.updateMol()

    def updateMol(self):
        self.mol = self.molecules[self.mlist.currentRow()]
        self.visual.setMol(self.mol)
        self.mtable.setMol(self.mol)
        self.mcell.setMol(self.mol)
        self.mkpoints.setMol(self.mol)
        self.pseglob.setMol(self.mol)
        self.psemol.setMol(self.mol)
        for tool in self.tools:
            tool.setMol(self.mol)
        self.fmt.blockSignals(True)
        self.fmt.setCurrentIndex(fmts.index(self.mol.getFmt()))
        self.fmt.blockSignals(False)
        undo = self.mol.getUndo()
        if undo:
            self.undoAction.setText('Undo '+undo)
            self.undoAction.setEnabled(True)
        else:
            self.undoAction.setText('Undo')
            self.undoAction.setDisabled(True)

    def updateParam(self):
        param = self.parameters[self.plist.currentRow()]
        self.param = param
        self.pedit.setParam(param)
