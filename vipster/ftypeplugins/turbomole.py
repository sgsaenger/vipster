# -*- coding: utf-8 -*-
from ..molecule import Molecule

name = 'turbomole'
extension = ''
argument = 'turbomole'

param = None

def parser(name,data):
    """ Parse turbomole specific input file """
    tmol = Molecule(name)
    tmol.setVec([[1,0,0],[0,1,0],[0,0,1]])
    tmol.setCellDim(40.0)
    for l in data[1:]:
        if '$' in l: break
        line = l.split()
        tmol.newAtom(line[3].upper(),line[0:3],'bohr')
    return tmol

writer = None
