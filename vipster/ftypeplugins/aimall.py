# -*- coding: utf-8 -*-
from ..molecule import Molecule

name = 'AIMALL'
extension = 'sum'
argument = 'aim'

param = None
writer = None

def parser(name,data):
    """ Parse AIMALL output to molecule

    Creates atoms from NACPs
    Other critical points will be parsed directly
    """
    tmol = Molecule(name)
    i=0
    while i<len(data):
        line = data[i].split()
        if line and line[0] == 'CP#':
            i+=1
            cptype = data[i].split()[3]
            if cptype == 'NACP':
                cptype = data[i].split()[4].rstrip('0123456789')
            tmol.newAtom(cptype,line[4:])
        i+=1
    return tmol,None
