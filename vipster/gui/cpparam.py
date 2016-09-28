# -*- coding: utf-8 -*-

from vipster.gui.qtwrapper import *


class CPParam(QWidget):
    def __init__(self):
        super(CPParam, self).__init__()
        self.updatedisable = False
        self.selector = QComboBox()
        self.selector.currentIndexChanged.connect(self.selectNameList)
        self.textArea = QTextEdit()
        self.textArea.textChanged.connect(self.saveNameList)
        vbox = QVBoxLayout()
        vbox.addWidget(self.selector)
        vbox.addWidget(self.textArea)
        vbox.setContentsMargins(0, 0, 0, 0)
        self.setLayout(vbox)

    def selectNameList(self):
        if self.updatedisable:
            return
        self.textArea.setText(self.param[str(self.selector.currentText())])

    def saveNameList(self):
        self.param[str(self.selector.currentText())] =\
            str(self.textArea.document().toPlainText())

    def setParam(self, param):
        self.updatedisable = True
        self.param = param
        self.selector.clear()
        self.selector.addItems([i for i in param if i[0] == '&'])
        self.updatedisable = False
