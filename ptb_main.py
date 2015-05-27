#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from sys import exit,argv,stdout
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
        backend.newMol()
    elif argv[1] == '-h':
        print_help(backend,0)
    elif argv[1] in backend.cli_indict:
        for i in argv[2:]:
            backend.readFile(argv[1],i,'cli')
    else:
        print_help(backend,1)

    app = QApplication([])
    gui = MakeWindow(backend,False)
    app.aboutToQuit.connect(gui.deleteLater)
    gui.loadView()
    exit(app.exec_())

def print_help(backend,err):
    f = stdout
    f.write('PWToolBox usage:\n')
    f.write('ptb_main [OPTIONS]\n\n')
    f.write('No option given: start GUI\n\n')
    f.write('Options:\n')
    f.write('-h: print this help\n')
    for i in backend.cli_indict.items():
        f.write(i[0]+' [FILES]: '+i[1].__doc__.split('\n')[0]+'\n')
    exit(err)

if __name__ == '__main__':
        main()
