# -*- coding: utf-8 -*-
from ..molecule import Molecule

name = 'PWScf Output'
extension = 'pwo'
argument = 'pwo'

param = None
writer = None

def parser(name,data):
    """ Parse PWScf output to trajectory """
    tmol = Molecule(name,steps=0)
    i=0
    vec=[[0,0,0],[0,0,0],[0,0,0]]
    gamma=False
    while i<len(data):
        line = data[i].split()
        #ignore empty lines
        if not line:
            pass
        #read number of atoms
        elif line[0:3] == ['number', 'of', 'atoms/cell']:
            nat = int(line[4])
        #read cell dimension
        elif line[0] == 'celldm(1)=':
            celldm = float(line[1])
        #read initial cell vectors
        elif line[0:2] == ['crystal','axes:']:
            for j in [0,1,2]:
                temp = data[i+1+j].split()
                vec[j]=[float(x) for x in temp[3:6]]
        #read initial positions:
        elif line[0] == 'site':
            tmol.newStep()
            tmol.set_celldm(celldm)
            tmol.set_vec(vec)
            for j in range(i+1,i+nat+1):
                atom = data[j].split()
                tmol.create_atom(atom[1],atom[6:9],'alat')
            i+=nat
        #read k-points:
        elif line[0] == 'gamma-point':
            gamma=True
        elif line[0:3] == ['number','of','k'] and not gamma:
            nk = int(line[4])
            kpoints=[]
            for j in range(i+2,i+nk+2):
                kp = data[j].split()
                kpoints.append([kp[4],kp[5],kp[6].strip('),'),kp[9]])
            tmol.set_kpoints('tpiba',kpoints)
            tmol.set_kpoints('active','tpiba')
            i+=nk
        #read step-vectors if cell is variable
        elif line[0] == 'CELL_PARAMETERS':
            for j in [0,1,2]:
                temp = data[i+1+j].split()
                vec[j]=[float(x) for x in temp[0:3]]
        #read step-coordinates
        elif line[0] == 'ATOMIC_POSITIONS':
            tmol.newStep()
            tmol.set_celldm(celldm)
            tmol.set_vec(vec)
            for j in range(i+1,i+nat+1):
                atom = data[j].split()
                tmol.create_atom(atom[0],atom[1:4],line[1].strip('()'),[int(x) for x in atom[4:]])
            i+=nat
        #break on reaching final coordinates (duplicate)
        elif line[0] == 'Begin':
            break
        #ignore everything else
        else:
            pass
        i+=1
    return tmol,None
