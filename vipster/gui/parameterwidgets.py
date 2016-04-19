# -*- coding: utf-8 -*-

from vipster.gui.cpparam import CPParam
from vipster.gui.pwparam import PWParam
from vipster.gui.lammpsparam import LammpsParam
from vipster.iowrapper import availParam,newParam

from PyQt4.QtGui import *

class ParamTab(QStackedWidget):
    def __init__(self):
        super(ParamTab,self).__init__()
        self.param = None
        self.paramList = [None,"pwi","cpmd","lmp"]
        self.addWidget(QLabel("No Parameter set selected"))
        self.addWidget(PWParam())
        self.addWidget(CPParam())
        self.addWidget(LammpsParam())

    def setParam(self,param):
        if param is None:
            self.setCurrentIndex(0)
        else:
            self.setCurrentIndex(self.paramList.index(param["type"]))
            self.currentWidget().setParam(param)

class ParamDialog(QDialog):
    def __init__(self,fmt,params):
        super(ParamDialog,self).__init__()
        self.label=QLabel("Please choose a valid parameter set:")
        self.selector=QComboBox()
        self.editor=ParamTab()
        self.selector.currentIndexChanged.connect(self.setParam)
        #setup parameter list:
        self.parameters = [newParam(fmt,i) for i in availParam(fmt)]
        self.parameters.extend([i for i in params if i["type"]==fmt])
        self.selector.addItems([i['name'] for i in self.parameters])
        #layout
        ok=QPushButton("OK")
        ok.clicked.connect(self.accept)
        cancel=QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        vbox=QVBoxLayout()
        vbox.addWidget(self.label)
        vbox.addWidget(self.selector)
        vbox.addWidget(self.editor)
        hbox=QHBoxLayout()
        hbox.addWidget(cancel)
        hbox.addWidget(ok)
        vbox.addLayout(hbox)
        self.setLayout(vbox)

    def setParam(self,idx):
        self.editor.setParam(self.parameters[self.selector.currentIndex()])

    def getParam(self):
        return self.parameters[self.selector.currentIndex()]
