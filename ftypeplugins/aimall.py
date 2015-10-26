# -*- coding: utf-8 -*-
from ptb.molecule import Trajectory

name = 'AIMALL'
extension = 'sum'
argument = 'aim'

writer = None

def parser(data):
    """ Parse AIMALL output to molecule

    Creates atoms from NACPs
    Other critical points will be parsed directly
    """
    tmol = Trajectory()
    i=0
    while i<len(data):
        line = data[i].split()
        if line and line[0] == 'CP#':
            i+=1
            cptype = data[i].split()[3]
            if cptype == 'NACP':
                cptype = data[i].split()[4].rstrip('0123456789')
            tmol.create_atom(cptype,line[4:],'bohr')
        i+=1
    return tmol,None
