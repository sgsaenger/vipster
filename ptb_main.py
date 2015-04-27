#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from sys import exit,argv
import cProfile
from ptb_mol import TBController
from ptb_gui import MakeWindow
from PyQt4.QtGui import QApplication

def main():
    """Main function

    Creates Backend, then parses arguments.
    If arguments match cli-interface, files will be parsed.
    If no cli-only behaviour is requested, gui will be started.
    """
    backend = TBController()

    if len(argv)==1:
        pass
    elif argv[1] == '-h':
        backend.print_help(0)
    elif argv[1] in backend.cli_indict:
        for i in argv[2:]:
            backend.readFile(argv[1],i,'cli')
    else:
        backend.print_help(1)

    app = QApplication([])
    gui = MakeWindow(backend,False)
    app.aboutToQuit.connect(gui.deleteLater)
    gui.loadView()
    exit(app.exec_())

if __name__ == '__main__':
        main()
