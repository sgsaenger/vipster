# -*- coding: utf-8 -*-
from ..molecule import Molecule

name = 'Empire XYZ'
extension = 'xyz'
argument = 'emp'

param = None

def parser(name,data):
    """ Parse Empire specific xyz file """
    tmol = Molecule(name)
    nat=int(data[0])
    tmol.setComment(data[1].strip())
    for j in range(2,nat+2):
        line=data[j].split()
        tmol.newAtom(line[0],line[1:4],'angstrom')
    vec=[0,0,0]
    for i in range(3):
        vec[i]=[float(x)for x in data[nat+3+i].split()]
    tmol.setVec(vec)
    tmol.setCellDim(1,fmt='angstrom')
    return tmol,None

def writer(mol,f,param,coordfmt):
    """
    Save Empire input file

    Same as xyz, cell parameters given after coordinates
    """
    # fixed format nat and comment
    f.write(str(mol.nat)+'\n')
    f.write('Hamil=PM3 calc=spt Periodic\n')
    # write coordinates
    for cntat in range(0,mol.nat):
            atom=mol.getAtom(cntat,'angstrom')
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                         atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
    f.write('\n')
    vec = mol.getVec()*mol.getCellDim('angstrom')
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[0]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[1]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[2]))
