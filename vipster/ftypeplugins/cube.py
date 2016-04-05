# -*- coding: utf-8 -*-
from .. import pse as glob_pse
from ..molecule import Molecule

name = 'Gaussian Cube'
extension = 'cub'
argument = 'cube'

param = None

def parser(name,data):
    """ Parse Gaussian Cube file """
    tmol = Molecule(name)
    fmt = "bohr"
    #parse data
    i=0
    #two lines of comments, combine
    tmol.comment = data[0]+";"+data[1]
    #nat, origin[3]
    nat=int(data[2].split()[0])
    origin=[float(data[2].split()[i]) for i in [1,2,3]]
    #n of datapoints(dir),cell_vec[3]
    nvol=[0,0,0]
    tvec=[0,0,0]
    for i in [0,1,2]:
        line=data[i+3].split()
        nvol[i]=int(line[0])
        if nvol[i]<0:
            nvol[i] = -nvol[i]
            fmt = "angstrom"
        tvec[i]=[float(line[j])*nvol[i] for j in [1,2,3]]
    tmol.setVec(tvec)
    tmol.nvol=nvol
    pse = list(glob_pse.keys())
    tmol.newAtoms(nat)
    for i in range(nat):
        # line = Z, charge(ignored), coord(x,y,z)
        line=data[i+6].split()
        tmol.setAtom(i,pse[int(line[0])],line[2:5])
    tmol.setFmt(fmt,scale=True)
    #rest of file has datagrid, x is outer loop, z inner
    tmol.setVol(nvol,data[6+nat:],origin)
    return tmol,None

writer = None
