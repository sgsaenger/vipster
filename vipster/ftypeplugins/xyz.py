# -*- coding: utf-8 -*-
from ..molecule import Molecule

name = 'xyz'
extension = 'xyz'
argument = 'xyz'
param = None

def parser(name,data):
    """ Parse xyz files (angstrom) """
    # create list of mol, trajectory support
    tmol = Molecule(name,steps=0)
    i=0
    while i < len(data):
            # handle empty lines at eof or between molecules
            if not data[i].strip().isdigit():
                    i+=1
                    continue
            #fixed format nat and comment
            tmol.newStep()
            nat = int(data[i])
            tmol.set_comment(data[i+1].strip())
            #read coordinates and types
            for j in range(i+2,i+nat+2):
                    line = data[j].split()
                    tmol.create_atom(line[0],line[1:4],'angstrom')
            i+=nat+2
    return tmol,None

def writer(mol,f,param,coordfmt):
    """ Save xyz in angstrom """
    # fixed format nat and comment
    f.write(str(mol.get_nat())+'\n')
    f.write(mol.get_comment()+'\n')
    # write coordinates
    for cntat in range(0,mol.get_nat()):
            atom=mol.get_atom(cntat,'angstrom')
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                         atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
    f.close()
