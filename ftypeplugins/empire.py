# -*- coding: utf-8 -*-
from ptb.molecule import Trajectory

name = 'Empire XYZ'
extension = 'xyz'
argument = 'emp'

param = None

def parser(name,data):
    """ Parse Empire specific xyz file """
    tmol = Trajectory(name,steps=1)
    nat=int(data[0])
    tmol.set_comment(data[1].strip())
    for j in range(2,nat+2):
        line=data[j].split()
        tmol.create_atom(line[0],line[1:4],'angstrom')
    vec=[0,0,0]
    for i in range(3):
        vec[i]=[float(x)for x in data[nat+3+i].split()]
    tmol.set_vec(vec)
    tmol.set_celldm(1,fmt='angstrom')
    return tmol,None

def writer(mol,f,param,coordfmt):
    """
    Save Empire input file

    Same as xyz, cell parameters given after coordinates
    """
    # fixed format nat and comment
    f.write(str(mol.get_nat())+'\n')
    f.write('Hamil=PM3 calc=spt Periodic\n')
    # write coordinates
    for cntat in range(0,mol.get_nat()):
            atom=mol.get_atom(cntat,'angstrom')
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                         atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
    f.write('\n')
    vec = mol.get_vec()*mol.get_celldm('angstrom')
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[0]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[1]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[2]))
