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
            tmol.comment = data[i+1].strip()
            #read coordinates and types
            tmol.newAtoms(nat)
            for j in range(nat):
                line=data[j+i+2].split()
                tmol.setAtom(j,line[0],line[1:4])
            tmol.setFmt('angstrom',scale=True)
            i+=nat+2
    return tmol,None

def writer(mol,f,param):
    """ Save xyz in angstrom """
    # fixed format nat and comment
    f.write(str(mol.nat)+'\n')
    f.write(mol.comment+'\n')
    # write coordinates
    for at in mol.getAtoms(fmt='angstrom'):
        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}\n'.format(at[0],*at[1]))
    f.close()
