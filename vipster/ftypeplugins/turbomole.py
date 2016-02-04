# -*- coding: utf-8 -*-
from ..molecule import Molecule

name = 'turbomole'
extension = 'turbo'
argument = 'turbo'

param = None

def parser(name,data):
    """ Parse turbomole specific input file """
    tmol = Molecule(name)
    for i in range(len(data)):
        if "$coord" in data[i]:
            for l in data[i+1:]:
                if '$' in l: break
                line = l.split()
                tmol.newAtom(line[3].capitalize(),line[0:3],'bohr')
            break
    return tmol,None

writer = None
