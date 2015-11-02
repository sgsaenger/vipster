# -*- coding: utf-8 -*-
import ptb
from ptb.molecule import Trajectory

name = 'Gaussian Cube file'
extension = 'cub'
argument = 'cube'

def parser(data):
    """ Parse Gaussian Cube file """
    tmol = Trajectory(steps=1)
    #parse data
    i=0
    #two lines of comments, combine
    tmol.set_comment(data[0]+";"+data[1])
    #nat, origin[3]
    nat=int(data[2].split()[0])
    origin=[float(data[2].split()[i]) for i in [1,2,3]]
    #n of datapoints(dir),cell_vec[3]
    nvol=[0,0,0]
    tvec=[0,0,0]
    for i in [0,1,2]:
        line=data[i+3].split()
        nvol[i]=int(line[0])
        tvec[i]=[float(line[j])*nvol[i] for j in [1,2,3]]
    tmol.set_vec(tvec)
    tmol.nvol=nvol
    pse = list(ptb.pse.keys())
    for i in range(nat):
        # line = Z, charge(ignored), coord(x,y,z)
        line=data[i+6].split()
        tmol.create_atom(pse[int(line[0])],line[2:5],'bohr')
    #rest of file has datagrid, x is outer loop, z inner
    tmol.set_vol(nvol,data[6+nat:],origin)
    return tmol,None

writer = None
#def writer(mol,f,param,coordfmt):
#    """ Save cube file """
#    f.write(mol.get_comment().split(';')[0])
#    f.write(mol.get_comment().split(';')[1])
#    vol = mol.get_vol()
#    s = vol.shape
#    vec = mol.get_vec()*mol.get_celldm()
#    vec = vec/s
#    f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(mol.get_nat(),0.,0.,0.))
#    f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(s[0],vec[0][0],vec[0][1],vec[0][2]))
#    f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(s[1],vec[1][0],vec[1][1],vec[1][2]))
#    f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(s[2],vec[2][0],vec[2][1],vec[2][2]))
#    for i in range(mol.get_nat()):
#        at = mol.get_atom(i,'bohr')
#        f.write('{:5.0f} {:.6f} {:.6f} {:.6f} {:.6f}\n'.format(mol.pse[at[0]][1],0,*at[1]))
#    vol = vol.flatten()
#    for i in range(1,vol.size+1):
#        f.write(str(vol[i-1])+'  ')
#        if i%6==0 or i%s[0]==0:
#            f.write('\n')
