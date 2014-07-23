#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from ptb_mol import TBController
from PyQt4.QtGui import *

#####################################################
# Application
#####################################################

def main():
        app = QApplication(sys.argv)
        control = TBController(sys.argv)
        sys.exit(app.exec_())

if __name__ == '__main__':
        main()
