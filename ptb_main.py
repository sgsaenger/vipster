#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import cProfile
from ptb_mol import TBController
from PyQt4.QtGui import *

#####################################################
# Application
#####################################################
def main():
        cProfile.run('main2()','stats.log')

def main2():
        app = TBController(sys.argv)
        sys.exit(app.exec_())

if __name__ == '__main__':
        main()
