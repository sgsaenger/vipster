# -*- coding: utf-8 -*-
from vipster.gui.qtwrapper import *

##############################################################
# Widget with titlebutton and collapsible main-widget
#############################################################


class collapsibleWidget(QWidget):
    def __init__(self, title, widget, show=True):
        super(collapsibleWidget, self).__init__()
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        button = QPushButton(title)
        button.clicked.connect(self.toggle)
        button.setFlat(True)
        self.widget = widget
        vlist = QVBoxLayout()
        vlist.addWidget(sep)
        vlist.addWidget(button)
        vlist.addWidget(self.widget)
        vlist.setContentsMargins(0, 0, 0, 0)
        self.setLayout(vlist)
        if not show:
            self.widget.hide()

    def toggle(self):
        if self.widget.isVisible():
            self.widget.hide()
        else:
            self.widget.show()
