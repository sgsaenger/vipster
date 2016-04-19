# -*- coding: utf-8 -*-
from vipster.molecule import Molecule

name = 'Empire XYZ'
extension = 'xyz'
argument = 'emp'

param = None

def parser(name,data):
    """ Parse Empire specific xyz file """
    tmol = Molecule(name)
    nat=int(data[0])
    tmol.comment = data[1].strip()
    tmol.newAtoms(nat)
    for j in range(nat):
        line=data[j+2].split()
        tmol.setAtom(j,line[0],line[1:4])
    tmol.setFmt('angstrom',scale=True)
    vec=[0,0,0]
    for i in range(3):
        vec[i]=[float(x)for x in data[nat+3+i].split()]
    tmol.setVec(vec)
    tmol.setCellDim(1,fmt='angstrom')
    return tmol,None

def writer(mol,f,param):
    """
    Save Empire input file

    Same as xyz, cell parameters given after coordinates
    """
    # fixed format nat and comment
    f.write(str(mol.nat)+'\n')
    f.write('Hamil=PM3 calc=spt Periodic\n')
    # write coordinates
    for at in mol.getAtoms(fmt='angstrom'):
        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}\n'.format(at[0],*at[1]))
    f.write('\n')
    vec = mol.getVec()*mol.getCellDim(fmt='angstrom')
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[0]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[1]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[2]))
