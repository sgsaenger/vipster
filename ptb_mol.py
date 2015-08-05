#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from math import sqrt
from collections import OrderedDict
from os.path import dirname,expanduser
from json import JSONDecoder

from molecule import Molecule,Trajectory
from ftypePlugins import cli_indict,gui_indict,gui_outdict

######################################################################
# MAIN CONTROLLER CLASS
######################################################################
class TBController(object):
    """
    I/O routines and handling of molecules/trajectories

    Dictionaries referencing read/write-routines
    Mol/Trajec and Param data saved in lists
    """

    def __init__(self):
            self._mol = []
            self._pwdata = []
            self._config = dict()
            self.pse = OrderedDict()
            self.cli_indict = cli_indict
            self.gui_indict = gui_indict
            self.gui_outdict= gui_outdict
            #self.cli_indict = OrderedDict([('-xyz',self._parseXyz),
            #                ('-pwi',self._parsePwi),
            #                ('-pwo',self._parsePwo),
            #                ('-pwof',self._parsePwoFinal),
            #                ('-lmp',self._parseLmp),
            #                ('-dmp',self._parseDmp),
            #                ('-cube',self._parseCube)])
            #self.indict = OrderedDict([('xyz',self._parseXyz),
            #               ('PWScf Input',self._parsePwi),
            #               ('PWScf Output' , self._parsePwo),
            #               ('PWO Final Conf.',self._parsePwoFinal),
            #               ('Gaussian Cube File',self._parseCube),
            #               ('Lammps Data File',self._parseLmp),
            #               ('Lammps Custom Dump',self._parseDmp)])
            #self.outdict= OrderedDict([('PWScf Input',self._writePwi),
            #               ('xyz',self._writeXyz),
            #               ('Empire xyz',self._writeEmpire),
            #               ('Gaussian Cube File',self._writeCube),
            #               ('Lammps Data File',self._writeLmp)])
            self.readConfig()

#####################################################################
# GET FUNCTIONS
#####################################################################

    def getMol(self,index):
        """ Return a given step of a given molecule """
        return self._mol[index]

    def getNMol(self):
        """ Return the number of loaded molecules/trajectories """
        return len(self._mol)

    def getPw(self,index):
        """ Return a given PW parameter set """
        return self._pwdata[index]

    def getNPw(self):
        """ Return the number of loaded parameter sets """
        return len(self._pwdata)

#####################################################################
# CREATE FUNCTIONS
#####################################################################

    def newMol(self):
        """ Create a new Trajectory with one empty Molecule """
        self._mol.append(Trajectory(self,1))

    def newTrajectory(self):
        """ Create a new (empty) Trajectory """
        self._mol.append(Trajectory(self,0))

    def newPw(self):
        """ Create a new (empty) dict for PWScf parameters """
        self._pwdata.append(OrderedDict())

#####################################################################
# READ FUNCTIONS
#####################################################################

    def readConfig(self):
        """Read config and PSE from json-file"""
        try:
            with open(expanduser('~/.toolbox.json')) as f:
                conf = JSONDecoder(object_pairs_hook=OrderedDict).decode(f.read())
        except:
            with open(dirname(__file__)+'/default.json') as f:
                conf = JSONDecoder(object_pairs_hook=OrderedDict).decode(f.read())
        self.pse=conf['PSE']
        self._config=conf['General']

    def readFile(self,fmt,filename,mode='gui'):
        """
        Read and parse a given file

        fmt -> file format, needs to be in indict/cli_indict
        filename -> path to file
        mode -> decide which dictionary to use

        If fmt is in dict, the file will be parsed
        """
        with open(filename,'r') as data:
            data = data.readlines()
            if mode =='cli':
                self.cli_indict[fmt](self,data)
            else:
                self.gui_indict[fmt](self,data)

    def _parseLmp(self,data):
        """ Parse Lammps data file

        Preliminary implementation!
        Element parsed via comment in 'Masses' part
        Only orthogonal cells supported
        Assumes angstrom
        """
        tmol=Molecule(self)
        i=0
        tvec=[[0,0,0],[0,0,0],[0,0,0]]
        while i< len(data):
            line = data[i].strip()
            if 'atoms' in line:
                nat = int(line.split()[0])
            elif 'atom types' in line:
                types = [0]*int(line.split()[0])
            elif 'xlo xhi' in line:
                tvec[0][0] = float(line.split()[1])-float(line.split()[0])
            elif 'ylo yhi' in line:
                tvec[1][1] = float(line.split()[1])-float(line.split()[0])
            elif 'zlo zhi' in line:
                tvec[2][2] = float(line.split()[1])-float(line.split()[0])
                tmol.set_vec(tvec)
                tmol.set_celldm(1,fmt='angstrom')
            elif 'Masses' in line:
                for j in range(i+2,i+2+len(types)):
                    if '#' in data[j]:
                        types[int(data[j].split()[0])-1]=data[j].split('#')[1].strip()
                    else:
                        raise NotImplementedError('cannot assign elements via masses')
                i+=len(types)+1
            elif 'Atoms' in line:
                for j in range(i+2,i+2+nat):
                    at = data[j].strip().split()
                    tmol.create_atom(types[int(at[1])-1],map(float,at[-3:]),'angstrom')
            i+=1
        self._mol.append([tmol])

    def _parseDmp(self,data):
        """ Parse Lammps custom dump files

        Preliminary implementation!
        Needs to be 'custom' format with either:
        'id element xs ys yz'
        or
        'id element x y z'
        Only orthogonal cells for now
        Assumes angstrom
        """
        tlist=[]
        i=0
        while i<len(data):
            line = data[i]
            if 'ITEM' in line:
                if 'TIMESTEP' in line:
                    i+=2
                    tmol=Molecule(self)
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
                        tmol.create_atom(at[1],map(float,at[2:]),fmt=fmt)
                    i+=nat+1
                    tlist.append(tmol)
            else:
                i+=1
        self._mol.append(tlist)

    def _parseCube(self,data):
        """ Parse Gaussian Cube file """
        tmol = Molecule(self)
        tcoord=[]
        tvec=[[0,0,0],[0,0,0],[0,0,0]]
        #parse data
        i=0
        #two lines of comments, combine
        tmol._comment=data[0]+";"+data[1]
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
        for i in range(nat):
            # line = Z, charge(ignored), coord(x,y,z)
            line=data[i+6].split()
            tmol.create_atom(self.pse.keys()[int(line[0])],map(float,line[2:5]),'bohr')
        #rest of file has datagrid, x is outer loop, z inner
        tmol.set_vol(nvol,data[6+nat:],origin)
        #finished molecule will be appended to list
        self._mol.append([tmol])

#############################################################################
# WRITE FUNCTIONS
#############################################################################

    def writeFile(self,ftype,mol,filename,param="",coordfmt=""):
        """
        Write a file to disk

        ftype -> type of file, needs to be in outdict
        mol -> molecule to save
        filename -> target filename
        param -> parameter set to save
        coordfmt -> format in which to save coordinates (bohr/angstrom/crystal/alat)
        """
        with open(filename,'w') as f:
            self.gui_outdict[ftype](mol,f,param,coordfmt)

    def _writeXyz(self,mol,f,param,coordfmt):
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

    def _writeEmpire(self,mol,f,param,coordfmt):
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
        vec = mol.get_vec()*mol.get_celldm()*0.52917721092
        f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[0]))
        f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[1]))
        f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[2]))

    def _writeCube(self,mol,f,param,coordfmt):
        f.write(mol.get_comment().split(';')[0])
        f.write(mol.get_comment().split(';')[1])
        vol = mol.get_vol()
        s = vol.shape
        vec = mol.get_vec()*mol.get_celldm()
        vec = vec/s
        f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(mol.get_nat(),0.,0.,0.))
        f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(s[0],vec[0][0],vec[0][1],vec[0][2]))
        f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(s[1],vec[1][0],vec[1][1],vec[1][2]))
        f.write('{:5d} {:.6f} {:.6f} {:.6f}\n'.format(s[2],vec[2][0],vec[2][1],vec[2][2]))
        for i in range(mol.get_nat()):
            at = mol.get_atom(i,'bohr')
            f.write('{:5.0f} {:.6f} {:.6f} {:.6f} {:.6f}\n'.format(mol.pse[at[0]][1],0,*at[1]))
        vol = vol.flatten()
        for i in range(1,vol.size+1):
            f.write(str(vol[i-1])+'  ')
            if i%6==0 or i%s[0]==0:
                f.write('\n')

    def _writeLmp(self,mol,f,param,coordfmt):
        """
        Save Lammps input data

        Preliminary implementation!
        Cell parameters have to be in lower triangular form
        Non orthogonal cell not supported for now
        Assumes angstrom until dialog is established
        """
        f.write('\n'+str(mol.get_nat())+' atoms\n')
        f.write(str(mol.get_ntyp())+' atom types\n\n')
        #check if box is orthogonal:
        v=mol.get_vec()*mol.get_celldm(fmt='angstrom')
        if not v.diagonal(1).any() and not v.diagonal(2).any():
            if not v.diagonal(-1).any() and not v.diagonal(-2).any():
                f.write('{:.5f} {:.5f} xlo xhi\n'.format(0.0,v[0][0]))
                f.write('{:.5f} {:.5f} ylo yhi\n'.format(0.0,v[1][1]))
                f.write('{:.5f} {:.5f} zlo zhi\n\n'.format(0.0,v[2][2]))
            else:
                raise TypeError('Non-orthogonal cell not yet supported')
        else:
            raise TypeError('Cell not in suitable Lammps format')
        f.write('Masses\n\n')
        t=list(mol.get_types())
        for i,j in enumerate(t):
            f.write('{:d} {:2.4f} #{:s}\n'.format(i+1,mol.pse[j][2],j))
        f.write('\nAtoms\n\n')
        for i in range(mol.get_nat()):
            at=mol.get_atom(i,'angstrom')
            f.write(('{:d} {:d}'+' {:d}'*1+' {:15.10f} {:15.10f} {:15.10f}\n').format(
                i,t.index(at[0])+1,0,*at[1]))
        f.write('\n')
