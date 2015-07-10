#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from math import sqrt
from collections import OrderedDict
from os.path import dirname,expanduser
from json import JSONDecoder

from molecule import Molecule

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
            self.cli_indict = OrderedDict([('-xyz',self._parseXyz),
                            ('-pwi',self._parsePwi),
                            ('-pwo',self._parsePwo),
                            ('-pwof',self._parsePwoFinal),
                            ('-lmp',self._parseLmp),
                            ('-dmp',self._parseDmp),
                            ('-cube',self._parseCube)])
            self.indict = OrderedDict([('xyz',self._parseXyz),
                           ('PWScf Input',self._parsePwi),
                           ('PWScf Output' , self._parsePwo),
                           ('PWO Final Conf.',self._parsePwoFinal),
                           ('Gaussian Cube File',self._parseCube),
                           ('Lammps Data File',self._parseLmp),
                           ('Lammps Custom Dump',self._parseDmp)])
            self.outdict= OrderedDict([('PWScf Input',self._writePwi),
                           ('xyz',self._writeXyz),
                           ('Empire xyz',self._writeEmpire),
                           ('Gaussian Cube File',self._writeCube),
                           ('Lammps Data File',self._writeLmp)])
            self.readConfig()

#####################################################################
# GET FUNCTIONS
#####################################################################

    def get_mol(self,index,step):
        """ Return a given step of a given molecule """
        return self._mol[index][step]

    def get_lmol(self,index):
        """ Return the length of a given trajectory """
        return len(self._mol[index])

    def get_nmol(self):
        """ Return the number of loaded molecules/trajectories """
        return len(self._mol)

    def get_pw(self,index):
        """ Return a given PW parameter set """
        return self._pwdata[index]

    def get_npw(self):
        """ Return the number of loaded parameter sets """
        return len(self._pwdata)

#####################################################################
# NEW MOLECULE
#####################################################################

    def newMol(self):
        """ Create a new (empty) Molecule """
        self._mol.append([Molecule(self)])

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
                self.cli_indict[fmt](data)
            else:
                self.indict[fmt](data)

    def _parseXyz(self,data):
        """ Parse xyz files (angstrom) """
        # create list of mol, trajectory support
        tlist = []
        # if first line does not contain nat exit
        if not data[0].strip().isdigit():
                print('not a xyz file')
                return
        i=0
        while i < len(data):
                # handle empty lines at eof or between molecules
                if not data[i].strip().isdigit():
                        i+=1
                        continue
                # create new molecule
                tmol = Molecule(self)
                #fixed format nat and comment
                nat = int(data[i])
                tmol._comment = data[i+1].strip()
                #read coordinates and types
                for j in range(i+2,i+nat+2):
                        line = data[j].split()
                        tmol.create_atom(line[0],map(float,line[1:4]),'angstrom')
                i+=nat+2
                tlist.append(tmol)
        #append to list of molecules
        self._mol.append(tlist)

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

    def _parsePwi(self,data):
        """ Parse PWScf input files

        Namelists will be parsed and saved in PW parameter set
        Supported Cards:
        - ATOMIC_SPECIES
        - ATOMIC_POSITIONS
        - K_POINTS
        - CELL_PARAMETERS
        Not supported:
        - CONSTRAINTS
        - OCCUPATIONS
        - ATOMIC_FORCES (PWSCFv5)
        """
        # no need for list, only one molecule per file
        tmol = Molecule(self)
        tparam = OrderedDict()
        tcoord = []
        tvec = [[0,0,0],[0,0,0],[0,0,0]]
        #parse data and create tparam
        while data:
            header = data.pop(0).strip().split()
            # ignore empty lines
            if not header:
                pass
                #debug output. case not really needed.
            # parse namelists
            elif header[0][0] == '&':
                tnl = OrderedDict()
                # parse entries
                line = data.pop(0).strip().split(',')
                while line[0] != '/':
                    for j in range(len(line)):
                        if line[j]:
                            tnl[line[j].split('=')[0].strip()]=line[j].split('=')[1].strip()
                    line = data.pop(0).strip().split(',')
                tparam[header[0].lower()]=tnl
            # parse card
            elif header[0][0].isupper():
                # 7 types of cards, need hardcoding
                # 4 supported for now

                #ATOMIC_SPECIES:
                #Name   Weight  PP-file
                if header[0] == 'ATOMIC_SPECIES':
                    for i in range(int(tparam['&system']['ntyp'])):
                        line = data.pop(0).strip().split()
                        tmol.pse[line[0]][2] = float(line[1])
                        tmol.pse[line[0]][0] = line[2]

                #ATOMIC_POSITIONS fmt
                #Name   x   y   z
                #Saved in temp list for parsing
                elif header[0] == 'ATOMIC_POSITIONS':
                    fmt = header[1].strip('()')
                    for i in range(int(tparam['&system']['nat'])):
                        #support empty lines
                        temp = data.pop(0).strip().split()
                        while not temp:
                                temp = data.pop(0).strip().split()
                        tcoord.append(temp)

                #K_POINTS fmt
                elif header[0] == 'K_POINTS':
                    #Gamma point only
                    if header[1] == 'gamma':
                        tmol.set_kpoints('active','gamma')
                    #MPGrid:
                    #x y z offset
                    #passed as whole string for now
                    elif header[1] == 'automatic':
                        line = data.pop(0).strip().split()
                        tmol.set_kpoints('automatic',line)
                        tmol.set_kpoints('active','automatic')
                    #else:
                    #number of kpoints
                    #x y z weight
                    #passed as whole string for now
                    else:
                        nk = int(data.pop(0).strip().split()[0])
                        kpoints = []
                        for i in range(nk):
                            kpoints.append(data.pop(0).strip().split())
                        tmol.set_kpoints('disc',kpoints)
                        tmol.set_kpoints('active',header[1])

                #CELL_PARAMETERS tbd
                #only needed if ibrav=0
                #tbd changed between pw4 and pw5, ignored for now
                elif header[0] == 'CELL_PARAMETERS':
                    for i in [0,1,2]:
                        line = data.pop(0).strip().split()
                        tvec[i]=[float(x)for x in line]
                else:
                    pass

        #Identify cell parameter representation, parse it
        #set celldm from input
        tmol.set_celldm(tparam['&system']['celldm(1)'])
        if tparam['&system']['ibrav'] == '0':
            #check if CELL_PARAMETERS card has been read, if not present, throw error
            if tvec == [[0,0,0],[0,0,0],[0,0,0]]:
                    print('CELL_PARAMETERS missing')
            else:
                    tmol.set_vec(tvec)
        elif tparam['&system']['ibrav'] == '1':
            #simple cubic is the default for new molecules
            #do nothing
            pass
        elif tparam['&system']['ibrav'] == '2':
            #face centered cubic
            tmol.set_vec([[-0.5,0,0.5],[0,0.5,0.5],[-0.5,0.5,0]])
        elif tparam['&system']['ibrav'] == '3':
            #body centered cubic
            tmol.set_vec([[0.5,0.5,0.5],[-0.5,0.5,0.5],[-0.5,-0.5,0.5]])
        elif tparam['&system']['ibrav'] == '4':
            #hexagonal
            if not 'celldm(3)' in tparam['&system']:
                print('celldm(3) missing')
            else:
                ca = float(tparam['&system']['celldm(3)'])
                tmol.set_vec([[1,0,0],[-0.5,sqrt(3)*0.5,0],[0,0,ca]])
        elif tparam['&system']['ibrav'] == '5':
            #trigonal
            if not 'celldm(4)' in tparam['&system']:
                print('celldm(4) missing')
            else:
                c = float(tparam['&system']['celldm(4)'])
                tx=sqrt((1-c)/2)
                ty=sqrt((1-c)/6)
                tz=sqrt((1+2*c)/3)
                tmol.set_vec([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]])
        elif tparam['&system']['ibrav'] == '-5':
            #trigonal,alternative
            if not 'celldm(4)' in tparam['&system']:
                print('celldm(4) missing')
            else:
                c = float(tparam['&system']['celldm(4)'])
                tx=sqrt((1-c)/2)
                ty=sqrt((1-c)/6)
                tz=sqrt((1+2*c)/3)
                u=(tz-2*sqrt(2)*ty)/sqrt(3)
                v=(tz+sqrt(2)*ty)/sqrt(3)
                tmol.set_vec([[u,v,v],[v,u,v],[v,v,u]])
        elif tparam['&system']['ibrav'] == '6':
            #simple tetragonal
            if not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                tmol.set_vec([[1,0,0],[0,1,0],[0,0,ca]])
        elif tparam['&system']['ibrav'] == '7':
            #body centered tetragonal
            if not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                tmol.set_vec([[0.5,-0.5,ca*0.5],[0.5,0.5,ca*0.5],[-0.5,-0.5,ca*0.5]])
        elif tparam['&system']['ibrav'] == '8':
            #simple orthorhombic
            if not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                tmol.set_vec([[1,0,0],[0,ba,0],[0,0,ca]])
        elif tparam['&system']['ibrav'] == '9':
            #basis centered orthorhombic
            if not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                tmol.set_vec([[0.5,ba*0.5,0],[-0.5,ba*0.5,0],[0,0,ca]])
        elif tparam['&system']['ibrav'] == '10':
            #face centered orthorhombic
            if not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                tmol.set_vec([[0.5,0,ca*0.5],[0.5,ba*0.5,0],[0,ba*0.5,ca*0.5]])
        elif tparam['&system']['ibrav'] == '11':
            #body centered orthorhombic
            if not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                tmol.set_vec([[0.5,ba*0.5,ca*0.5],[-0.5,ba*0.5,ca*0.5],[-0.5,-ba*0.5,ca*0.5]])
        elif tparam['&system']['ibrav'] == '12':
            #simple monoclinic
            if not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            elif not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(4)' in tparam['&system']:
                print(  'celldm(4) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                cg = float(tparam['&system']['celldm(4)'])
                tmol.set_vec([[1,0,0],[ba*cg,ba*sqrt(1-cg),0],[0,0,ca]])
        elif tparam['&system']['ibrav'] == '-12':
            #simple monoclinic, alternate definition
            if not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            elif not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(5)' in tparam['&system']:
                print(  'celldm(5) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                cb = float(tparam['&system']['celldm(5)'])
                tmol.set_vec([[1,0,0],[0,ba,0],[ca*cb,0,ca*sqrt(1-cb)]])
        elif tparam['&system']['ibrav'] == '13':
            #base centered monoclinic
            if not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            elif not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(4)' in tparam['&system']:
                print(  'celldm(4) missing' )
            else:
                ca = float(tparam['&system']['celldm(3)'])
                ba = float(tparam['&system']['celldm(2)'])
                cg = float(tparam['&system']['celldm(4)'])
                tmol.set_vec([[0.5,0,-ca*0.5],[ba*cg,ba*sqrt(1-cg),0],[0.5,0,ca*0.5]])
        elif tparam['&system']['ibrav'] == '14':
            #base centered monoclinic
            if not 'celldm(2)' in tparam['&system']:
                print(  'celldm(2) missing' )
            elif not 'celldm(3)' in tparam['&system']:
                print(  'celldm(3) missing' )
            elif not 'celldm(4)' in tparam['&system']:
                print(  'celldm(4) missing' )
            elif not 'celldm(5)' in tparam['&system']:
                print(  'celldm(5) missing' )
            elif not 'celldm(6)' in tparam['&system']:
                print(  'celldm(6) missing' )
            else:
                ba = float(tparam['&system']['celldm(2)'])
                ca = float(tparam['&system']['celldm(3)'])
                cg = float(tparam['&system']['celldm(4)'])
                cb = float(tparam['&system']['celldm(5)'])
                cal = float(tparam['&system']['celldm(6)'])
                tmol.set_vec([[1,0,0],[ba*cg,ba*sqrt(1-cg),0],
                        [ca*cb,ca*(cal-cb*cg)/sqrt(1-cg),ca*sqrt(1+2*cal*cb*cg-cal*cal-cb*cb-cg*cg)/sqrt(1-cg)]])


        #create atoms after creating cell:
        for i in range(len(tcoord)):
            if len(tcoord[i])>4:
                tmol.create_atom(tcoord[i][0],map(float,tcoord[i][1:4]),fmt,
                        map(int,tcoord[i][4:]))
            else:
                tmol.create_atom(tcoord[i][0],map(float,tcoord[i][1:4]),fmt)
        #delete nat, ntype and celldm before returning to controller
        del tparam['&system']['nat']
        del tparam['&system']['ntyp']
        for i in range(1,7):
            test='celldm('+str(i)+')'
            if test in tparam['&system']:
                del tparam['&system'][test]

        #Append to controller
        self._mol.append([tmol])
        self._pwdata.append(tparam)

    def _parsePwo(self,data):
        """ Parse PWScf output to trajectory """
        #Multiple configs supported
        tlist = []
        #read list of molecules:
        i=0
        vec=[[0,0,0],[0,0,0],[0,0,0]]
        while i<len(data):
            line = data[i].split()
            #ignore empty lines
            if not line:
                pass
            #read number of atoms
            elif line[0:3] == ['number', 'of', 'atoms/cell']:
                nat = int(line[4])
            #TODO: tweak to recognize celldm(n), save if !0 in param
            #read cell dimension
            elif line[0] == 'celldm(1)=':
                celldm = float(line[1])
            #TODO: care for different formats
            #read initial cell vectors
            elif line[0:2] == ['crystal','axes:']:
                for j in [0,1,2]:
                    temp = data[i+1+j].split()
                    vec[j]=[float(x) for x in temp[3:6]]
            # read initial positions:
            elif line[0] == 'site':
                tmol = Molecule(self)
                tmol.set_celldm(celldm)
                tmol.set_vec(vec)
                for j in range(i+1,i+nat+1):
                    atom = data[j].split()
                    tmol.create_atom(atom[1],map(float,atom[6:9]),'alat')
                i+=nat
                tlist.append(tmol)
            #read step-vectors if cell is variable
            elif line[0] == 'CELL_PARAMETERS':
                for j in [0,1,2]:
                    temp = data[i+1+j].split()
                    vec[j]=[float(x) for x in temp[0:3]]
            #read step-coordinates
            elif line[0] == 'ATOMIC_POSITIONS':
                tmol = Molecule(self)
                tmol.set_celldm(celldm)
                tmol.set_vec(vec)
                for j in range(i+1,i+nat+1):
                    atom = data[j].split()
                    if len(atom)>4:
                        tmol.create_atom(atom[0],map(float,atom[1:4]),line[1].strip('()'),map(int,atom[4:]))
                    else:
                        tmol.create_atom(atom[0],map(float,atom[1:4]),line[1].strip('()'))
                i+=nat
                tlist.append(tmol)
            #break on reaching final coordinates (duplicate)
            elif line[0] == 'Begin':
                break
            #ignore everything else
            else:
                pass
            i+=1
        self._mol.append(tlist)

    def _parsePwoFinal(self,data):
        """ Parse only final configuration of PWScf output """
        i=0
        vec=[[0,0,0],[0,0,0],[0,0,0]]
        final = False
        while i<len(data):
            line = data[i].split()
            if not line:
                pass
            elif line[0:3] == ['number', 'of', 'atoms/cell']:
                nat = int(line[4])
            elif line[0] == 'celldm(1)=':
                celldm = float(line[1])
            elif line[0:2] == ['crystal','axes:']:
                for j in [0,1,2]:
                    temp = data[i+1+j].split()
                    vec[j]=[float(x) for x in temp[3:6]]
            elif line[0] == 'Begin':
                final = True
            #update cell parameters for vc calculations
            elif final and line[0] == 'CELL_PARAMETERS':
                for j in [0,1,2]:
                    temp = data[i+1+j].split()
                    vec[j]=[float(x) for x in temp[0:3]]

            elif final and line[0] == 'ATOMIC_POSITIONS':
                tmol = Molecule(self)
                tmol.set_celldm(celldm)
                tmol.set_vec(vec)
                for j in range(i+1,i+nat+1):
                    atom = data[j].split()
                    if len(atom)>4:
                        tmol.create_atom(atom[0],map(float,atom[1:4]),line[1].strip('()'),map(int,atom[4:]))
                    else:
                        tmol.create_atom(atom[0],map(float,atom[1:4]),line[1].strip('()'))
            else:
                pass
            i+=1
        self._mol.append([tmol])

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
        tmol.set_vol_gradient()
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
            self.outdict[ftype](mol,f,param,coordfmt)

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
        Format taken from coordfmt
        """
        f.write('\n'+str(mol.get_nat())+' atoms\n')
        f.write(str(mol.get_ntyp())+' atom types\n\n')
        #check if box is orthogonal:
        v=mol.get_vec()*mol.get_celldm(fmt=coordfmt)
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
            at=mol.get_atom(i,coordfmt)
            f.write(('{:d} {:d}'+' {:d}'*1+' {:15.10f} {:15.10f} {:15.10f}\n').format(
                i,t.index(at[0])+1,0,*at[1]))
        f.write('\n')

    def _writePwi(self,mol,f,param,coordfmt):
        """
        Save PWScf input file

        Needs both mol and param
        Respects coordfmt
        """
        #&control, &system and &electron namelists are mandatory
        for i in ['&control','&system','&electrons']:
            f.write(i+'\n')
            #write all parameters
            if i == '&system':
                f.write(' nat='+str(mol.get_nat())+'\n')
                f.write(' ntyp='+str(mol.get_ntyp())+'\n')
                f.write(' celldm(1)='+str(mol.get_celldm())+'\n')
            for j in range(len(param[i])):
                f.write(' '+param[i].keys()[j]+'='+param[i].values()[j]+'\n')
            f.write('/\n\n')
        #&ions only when needed
        if param['&control']['calculation'] in ["'relax'","'vc-relax'","'md'","'vc-md'"]:
            f.write('&ions'+'\n')
            for j in range(len(param['&ions'])):
                f.write(' '+param['&ions'].keys()[j]+'='+param['&ions'].values()[j]+'\n')
            f.write('/\n\n')
        #&cell only when needed
        if param['&control']['calculation'] in ["'vc-relax'","'vc-md'"]:
            f.write('&cell'+'\n')
            for j in range(len(param['&cell'])):
                f.write(' '+param['&cell'].keys()[j]+'='+param['&cell'].values()[j]+'\n')
            f.write('/\n\n')

        #ATOMIC_SPECIES card:
        f.write('ATOMIC_SPECIES'+'\n')
        types = list(mol.get_types())
        for i in range(len(mol.get_types())):
            atom = types[i]
            f.write(atom+'    '+str(mol.pse[atom][2])+'   '+str(mol.pse[atom][0])+'\n')
        f.write('\n')

        #ATOMIC_POSITIONS
        f.write('ATOMIC_POSITIONS'+' '+coordfmt+'\n')
        for i in range(mol.get_nat()):
            atom=mol.get_atom(i,coordfmt)
            if 0 in atom[3]:
                f.write('{:4s} {:15.10f} {:15.10f} {:15.10f} {:1d} {:1d} {:1d}'.format(
                    atom[0],atom[1][0],atom[1][1],atom[1][2],atom[3][0],atom[3][1],atom[3][2])+'\n')
            else:
                f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                    atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
        f.write('\n')

        #K_POINTS
        f.write('K_POINTS'+' '+mol.get_kpoints('active')+'\n')
        #Gamma point only
        if mol.get_kpoints('active') == 'gamma':
            pass
        #MPGrid:
        #x y z offset
        elif mol.get_kpoints('active') == 'automatic':
            auto = mol.get_kpoints('automatic')
            f.write('{:4s}{:4s}{:4s}{:4s}{:4s}{:4s}'.format(
                    auto[0],auto[1],auto[2],auto[3],auto[4],auto[5])+'\n')
        #number of kpoints
        #x y z weight
        else:
            disc=mol.get_kpoints('disc')
            f.write(str(len(disc))+'\n')
            for i in range(len(disc)):
                f.write('{:4s}{:4s}{:4s}{:4s}'.format(
                        disc[i][0],disc[i][1],disc[i][2],disc[i][3])+'\n')
        f.write('\n')

        #Cell parameters
        f.write('CELL_PARAMETERS'+'\n')
        fmt='{0[0][0]:15.10f} {0[0][1]:15.10f} {0[0][2]:15.10f}\n' + \
            '{0[1][0]:15.10f} {0[1][1]:15.10f} {0[1][2]:15.10f}\n' + \
            '{0[2][0]:15.10f} {0[2][1]:15.10f} {0[2][2]:15.10f}\n'
        f.write(fmt.format(mol.get_vec()))
