# -*- coding: utf-8 -*-
from vipster.molecule import Molecule

name = "Lammps Custom Dump"
argument = 'dmp'
extension = 'dmp'

param = None
writer = None


def parser(name, data):
    """ Parse Lammps custom dump files

    Preliminary implementation!
    Needs to be 'custom' format with either:
    'id element xs ys yz'
    or
    'id element x y z'
    Only orthogonal cells for now
    Assumes angstrom
    """
    tmol = Molecule(name, steps=0)
    i = 0
    while i < len(data):
        line = data[i]
        if 'ITEM' in line:
            if 'TIMESTEP' in line:
                i += 2
                tmol.newStep()
            elif 'NUMBER OF ATOMS' in line:
                nat = int(data[i + 1])
                i += 2
            elif 'BOX BOUNDS' in line:
                tvec = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
                tvec[0][0] = float(data[i + 1].split()[1]) -\
                    float(data[i + 1].split()[0])
                tvec[1][1] = float(data[i + 2].split()[1]) -\
                    float(data[i + 2].split()[0])
                tvec[2][2] = float(data[i + 3].split()[1]) -\
                    float(data[i + 3].split()[0])
                tmol.setVec(tvec)
                tmol.setCellDim(1, fmt='angstrom')
                i += 4
            elif 'ATOMS' in line:
                line = line.split()
                if 'id' not in line or 'element' not in line or\
                        ('xs' not in line and 'x' not in line):
                    raise NotImplementedError("Lammps dump in not (yet) "
                                              "recognized format")
                # ididx = line.index('id') - 2
                elidx = line.index('element') - 2
                if 'xs' in line:
                    xidx = line.index('xs') - 2
                    yidx = line.index('ys') - 2
                    zidx = line.index('zs') - 2
                    fmt = 'crystal'
                else:
                    xidx = line.index('x') - 2
                    yidx = line.index('y') - 2
                    zidx = line.index('z') - 2
                    fmt = 'angstrom'
                if 'q' in line:
                    qidx = line.index('q') - 2
                else:
                    qidx = False
                tmol.newAtoms(nat)
                for j in range(nat):
                    at = data[j + 1 + i].split()
                    if qidx:
                        tmol.setAtom(j, at[elidx],
                                     [at[xidx], at[yidx], at[zidx]],
                                     charge=float(at[qidx]))
                    else:
                        tmol.setAtom(j, at[elidx],
                                     [at[xidx], at[yidx], at[zidx]])
                tmol.setFmt(fmt, scale=True)
                i += nat + 1
        else:
            i += 1
    return tmol, None
