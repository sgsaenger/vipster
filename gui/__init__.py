# -*- coding: utf-8 -*-
from PyQt4.QtGui import QApplication
from .mainwidget import MainWidget

def launch(m=[],p=[]):
    app = QApplication([])
    gui = MainWidget(m,p)
    app.aboutToQuit.connect(gui.deleteLater)
    exit(app.exec_())
