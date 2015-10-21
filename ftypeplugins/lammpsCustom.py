# -*- coding: utf-8 -*-
name = "Lammps Custom Dump"
argument = '-dmp'
extension = 'dmp'

writer = None

def parser(controller,data):
    """ Parse Lammps custom dump files

    Preliminary implementation!
    Needs to be 'custom' format with either:
    'id element xs ys yz'
    or
    'id element x y z'
    Only orthogonal cells for now
    Assumes angstrom
    """
    controller.newTrajectory()
    tmol = controller.getMol(-1)
    i=0
    while i<len(data):
        line = data[i]
        if 'ITEM' in line:
            if 'TIMESTEP' in line:
                i+=2
                tmol.newMol()
            elif 'NUMBER OF ATOMS' in line:
                nat=int(data[i+1])
                i+=2
            elif 'BOX BOUNDS' in line:
                tvec=[[0,0,0],[0,0,0],[0,0,0]]
                tvec[0][0] = float(data[i+1].split()[1])-float(data[i+1].split()[0])
                tvec[1][1] = float(data[i+2].split()[1])-float(data[i+2].split()[0])
                tvec[2][2] = float(data[i+3].split()[1])-float(data[i+3].split()[0])
                tmol.set_vec(tvec)
                tmol.set_celldm(1,fmt='angstrom')
                i+=4
            elif 'ATOMS' in line:
                if line.strip() == 'ITEM: ATOMS id element xs ys zs':
                    fmt='crystal'
                elif line.strip() == 'ITEM: ATOMS id element x y z':
                    fmt='angstrom'
                else:
                    raise NotImplementedError('Lammps dump in not (yet) recognized format')
                for j in range(i+1,i+1+nat):
                    at = data[j].split()
                    tmol.create_atom(at[1],at[2:],fmt=fmt)
                i+=nat+1
        else:
            i+=1
