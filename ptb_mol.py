#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import copy
from ptb_gui import MainWindow
from math import sqrt,floor
from collections import OrderedDict
import numpy as np
from PyQt4.QtGui import *

######################################################################
# PSE DICTIONARY
######################################################################
# pse[0] id
# pse[1] ~ weight
pse={"X":  [0,0.0],

     "H" : [1,1.0079],
     "He": [2,4.0026],

     "Li": [3, 6.941],
     "Be": [4, 9.0122],
     "B" : [5,10.811],
     "C" : [6,12.0107],
     "N" : [7,14.007],
     "O" : [8,15.999],
     "F" : [9,18.998],
     "Ne": [10,20.18],

     "Na": [11,22.99],
     "Mg": [12,24.305],
     "Al": [13,26.982],
     "Si": [14,28.086],
     "P" : [15,30.974],
     "S" : [16,32.065],
     "Cl": [17,35.453],
     "Ar": [18,39.948],

     "K" : [19,39.098],
     "Ca": [20,40.078],
     "Ga": [31,69.723],
     "Ge": [32,72.64],
     "As": [33,74.922],
     "Se": [34,78.96],
     "Br": [35,79.904],
     "Kr": [36,83.798]

     }

######################################################################
# MAIN CONTROLLER CLASS
######################################################################
class TBController:

        def __init__(self,argv):
                self.mol = []
                self.pwdata = []
                #self.data = []
                self.gui = MainWindow(self)
                self.indict = OrderedDict([('xyz',self.parseXyz),
                               ('PWScf Input',self.parsePwi),
                               ('PWScf Output' , self.parsePwo)])
                self.outdict= OrderedDict([('PWScf Input',self.writePwi),
                               ('xyz',self.writeXyz)])
                self.argumentHandler(argv)


#####################################################################
# Handle command line arguments:
#####################################################################
        def argumentHandler(self,argv):
                for i in range(1,len(argv)):
                        if argv[i] == '-pwi':
                                self.readFile('PWScf Input',argv[i+1])
                                self.gui.centralWidget().loadView()
                        elif argv[i] == '-xyz':
                                self.readFile('xyz',argv[i+1])
                                self.gui.centralWidget().loadView()


#####################################################################
# GET FUNCTIONS
#####################################################################

        def get_mol(self,index,step):
                return self.mol[index][step]

        def get_lmol(self,index):
                return len(self.mol[index])

        def get_nmol(self):
                return len(self.mol)

        def get_pw(self,index):
                return self.pwdata[index]

        def get_npw(self):
                return len(self.pwdata)


#####################################################################
# READ FUNCTIONS
#####################################################################

        def readFile(self,fmt,filename):
                data = open(filename,'r').readlines()
                self.indict[fmt](data)
                #return data

        def parseXyz(self,data):
                # create list of mol, trajectory support
                tlist = []
                # if first line does not contain nat exit
                if not data[0].strip().isdigit():
                        print 'not a xyz file'
                        return
                i=0
                while i < len(data):
                        # handle empty lines at eof or between molecules
                        if not data[i].strip().isdigit(): i+=1
                        # create new molecule
                        tmol = Molecule()
                        #fixed format nat and comment
                        nat = int(data[i])
                        tmol.comment = data[i+1]
                        #read coordinates and types
                        for j in range(i+2,i+nat+2):
                                line = data[j].split()
                                tmol.create_atom(line[0],float(line[1]),float(line[2]),float(line[3]))
                        i+=nat+2
                        tlist.append(tmol)
                #append to list of molecules
                self.mol.append(tlist)

        def parsePwi(self,data):
                # no need for list, only one molecule per file
                tmol = Molecule()
                tparam = PWParam()
                tcoord = []
                #parse data and create tparam
                while data:
                        header = data.pop(0).strip().split()
                        # ignore empty lines
                        if not header:
                                pass
                                #debug output. case not really needed.
                        # parse namelists
                        elif header[0][0] == '&':
                                tnl = {}
                                # parse entries
                                line = data.pop(0).strip().split(',')
                                while line[0] != '/':
                                        for j in range(len(line)):
                                                tnl[line[j].split('=')[0].strip()]=line[j].split('=')[1].strip()
                                        line = data.pop(0).strip().split(',')
                                tparam[header[0]]=tnl
                                #debug output
                        # parse card
                        elif header[0][0].isupper():
                                # 7 types of cards, need hardcoding
                                # 4 supported for now

                                #ATOMIC_SPECIES:
                                #Name   Weight  PP-file
                                if header[0] == 'ATOMIC_SPECIES':
                                        for i in range(int(tparam['&system']['ntyp'])):
                                                line = data.pop(0).strip().split()
                                                tparam['pse'][line[0]][1] = float(line[1])
                                                tparam['pse'][line[0]].append(line[2])

                                #ATOMIC_POSITIONS fmt
                                #Name   x   y   z
                                #Saved in temp list for parsing
                                elif header[0] == 'ATOMIC_POSITIONS':
                                        fmt = header[1]
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
                                                tparam['K_POINTS']=['gamma']
                                        #MPGrid:
                                        #x y z offset
                                        #passed as whole string for now
                                        elif header[1] == 'automatic':
                                                #line = data.pop(0).strip().split()
                                                line = data.pop(0)
                                                #tparam['K_POINTS']=['automatic',[line[0:3],line[3:7]]]
                                                tparam['K_POINTS']=['automatic',line]
                                        #else:
                                        #number of kpoints
                                        #x y z weight
                                        #passed as whole string for now
                                        else:
                                                nk = data.pop(0).strip().split()
                                                kpoints = []
                                                for i in range(nk):
                                                        #kpoints.append(data.pop(0).strip().split())
                                                        kpoints.append(data.pop(0))
                                                taparam['K_POINTS']=[header[1],nk,kpoints]

                                elif header[0] == 'CELL_PARAMETERS':
                                        vec=[[0,0,0],[0,0,0],[0,0,0]]
                                        for i in range(3):
                                                line = data.pop(0).strip().split()
                                                vec[i]=[float(x)for x in line]
                                        tmol.set_vec(vec)
                                else:
                                        pass

                #Create Molecule after getting the cell parameters
                #set celldm from input
                tmol.set_celldm(tparam['&system']['celldm(1)'])
                #create atoms:
                for i in range(len(tcoord)):
                        tmol.create_atom(tcoord[i][0],float(tcoord[i][1]),float(tcoord[i][2]),float(tcoord[i][3]),fmt)
                #delete nat and ntype before returning to controller
                del tparam['&system']['nat']
                del tparam['&system']['ntyp']

                #Append to controller
                self.mol.append([tmol])
                self.pwdata.append(tparam)

        def parsePwo(self,data):
                #Multiple configs supported
                tlist = []
                #no parameters for now
                #tparam=PWParam()
                #read list of molecules:
                i=0
                while i<len(data):
                        line = data[i].split()
                        if not line: 
                                pass
                        elif line[0:3] == ['number', 'of', 'atoms/cell']:
                                nat = int(line[4])
                        elif line[0] == 'celldm(1)=':
                                celldm = float(line[1])
                        elif line[0:2] == ['crystal','axes:']:
                                vec=[[0,0,0],[0,0,0],[0,0,0]]
                                for j in range(3):
                                        temp = data[i+1+j].split()
                                        vec[j]=[float(x) for x in temp[3:6]]
                        elif line[0] == 'ATOMIC_POSITIONS':
                                tmol = Molecule()
                                tmol.set_celldm(celldm)
                                tmol.set_vec(vec)
                                for j in range(i+1,i+nat+1):
                                        atom = data[j].split()
                                        tmol.create_atom(atom[0],float(atom[1]),float(atom[2]),float(atom[3]),line[1].strip(')').strip('('))
                                i+=nat
                                tlist.append(tmol)
                        else:
                                pass
                        i+=1
                self.mol.append(tlist)


#############################################################################
# WRITE FUNCTIONS
#############################################################################

        def writeFile(self,ftype,mol,filename,param="",coordfmt=""):
                self.outdict[ftype](mol,filename,param,coordfmt)

        def writeXyz(self,mol,filename,param,coordfmt):
                if filename == "":
                        f=sys.stdout
                else:
                        f=open(filename,'w')
                # fixed format nat and comment
                f.write(str(mol.get_nat())+'\n')
                f.write(mol.get_comment()+'\n')
                # write coordinates
                for cntat in range(0,mol.get_nat()):
                        atom=mol.get_atom(cntat,'angstrom')
                        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                                     atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
                f.close()

        def writePwi(self,mol,filename,param,coordfmt):
                if filename == "":
                        f=sys.stdout
                else:
                        f=open(filename,'w')

                #&control, &system and &electron namelists are mandatory
                for i in ['&control','&system','&electrons']:
                        f.write(i+'\n')
                        #write all parameters
                        if i == '&system':
                                f.write(' nat='+str(mol.get_nat())+'\n')
                                f.write(' ntyp='+str(mol.get_ntyp())+'\n')
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
                        f.write(atom+'    '+str(param['pse'][atom][1])+'   '+str(param['pse'][atom][2])+'\n')
                f.write('\n')

                #ATOMIC_POSITIONS
                f.write('ATOMIC_POSITIONS'+' '+coordfmt+'\n')
                for i in range(mol.get_nat()):
                        atom=mol.get_atom(i)
                        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                            atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
                f.write('\n')

                #K_POINTS
                f.write('K_POINTS'+' '+param['K_POINTS'][0]+'\n')
                #Gamma point only
                if param['K_POINTS'][0] == 'gamma':
                        pass
                #MPGrid:
                #x y z offset
                #passed as whole string for now
                elif param['K_POINTS'][0] == 'automatic':
                        #f.write('{:4i}{:4i}{:4i}{:4i}{:4i}{:4i}'.format(
                        #        param['K_POINTS'][1][0][0],
                        #        param['K_POINTS'][1][0][1],
                        #        param['K_POINTS'][1][0][2],
                        #        param['K_POINTS'][1][0][0],
                        #        param['K_POINTS'][1][1][1],
                        #        param['K_POINTS'][1][1][2]
                        #        )+'\n\n')
                        f.write(param['K_POINTS'][1])
                #number of kpoints
                #x y z weight
                #passed as whole string for now
                #UNTESTED, beware!
                else:
                        f.write(param['K_POINTS'][1])
                        for i in range(param['K_POINTS'][1]):
                                f.write(param['K_POINTS'][2][i])
                f.write('\n')

                #Cell parameters
                f.write('CELL_PARAMETERS'+'\n')
                vec = mol.get_vec()
                fmt='{0[0][0]:15.10f} {0[0][1]:15.10f} {0[0][2]:15.10f}\n' + \
                    '{0[1][0]:15.10f} {0[1][1]:15.10f} {0[1][2]:15.10f}\n' + \
                    '{0[2][0]:15.10f} {0[2][1]:15.10f} {0[2][2]:15.10f}\n'
                f.write(fmt.format(mol.get_vec()))

                #Close file
                f.close()

######################################################################
# MOLECULE CLASS
######################################################################
class Molecule:

        def __init__(self):
                # set atom list
                self.at_n=[]
                self.at_c=[]
                self.celldm = 1.0
                self.vec=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self.vecinv=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self.comment = ''

        ######################################################
        # ATOM FUNCTIONS
        ######################################################

        # append new atom
        def create_atom(self,name='C',x=0.,y=0.,z=0.,fmt='angstrom'):
                #self.at.append(self.Atom(self,name,x,y,z,fmt))
                self.at_n.append(name)
                self.at_c.append(self.set_coord([x,y,z],fmt))

        # append copy of existing atom
        def append_atom_cp(self,addat):
                #self.at.append(copy.copy(addat))
                self.at_n.append(addat[0])
                self.at_c.append(self.set_coord(addat[1],addat[2]))

        # insert atom at given position
        def insert_atom(self,pos,addat):
                #self.at.insert(pos,copy.copy(addat))
                self.at_n.insert(pos,adddat[0])
                self.at_c.insert(pos,self.set_coord(addat[1],addat[2]))

        # remove atom
        def del_atom(self,index):
                del self.at[index]

        # append molecule
        #def append_mol(self, mol):
        #        for i in range(mol.get_nat()):
        #                self.at.append(copy.copy(mol.at[i]))

        ######################################################
        # CELL MULTIPLICATION
        ######################################################

        #TODO
        def getOffsets(self,mult):
                vec = self.get_vec()*self.celldm
                cent = self.get_center()
                off = []
                tmult = [1,1,1]
                #save the multiplicators for vec:
                for i in [0,1,2]:
                        if mult[i]%2 == 0:
                                tmult[i]=[x+0.5-mult[i]/2 for x in range(mult[i])]
                        else:
                                tmult[i]=[x-floor(mult[i]/2) for x in range(mult[i])]
                #generate offsets:
                for i in tmult[0]:
                        for j in tmult[1]:
                                for k in tmult[2]:
                                        off.append((i*vec[0]+j*vec[1]+k*vec[2])-cent)
                return off

        ######################################################
        # SET FUNCTIONS
        ######################################################

        def set_atom(self,index,name,coord,fmt):
                self.at_n[index]=name
                self.at_c[index]=self.set_coord(coord,fmt)

        def set_comment(self,comment):
                self.comment = comment

        # set celldm
        def set_celldm(self,cdm):
                self.celldm = float(cdm)
                self.set_center()

        # set vectors
        def set_vec(self,vec):
                self.vec = np.array(vec).T
                self.vecinv = np.linalg.inv(self.vec)
                self.set_center()

        # set center of cell
        def set_center(self):
                vec = self.get_vec()
                x = self.vec[0]
                y = self.vec[1]
                z = self.vec[2]
                self.center = (vec[0]+vec[1]+vec[2])*self.celldm/2

        #######################################################
        # COORD FMT FUNCTIONS
        # to be called only by atom set/get
        ######################################################

        def set_coord(self,coord,fmt='bohr'):
                coord = np.array(coord)
                if fmt == 'angstrom':
                        return coord/0.52917721092
                elif fmt == 'bohr':
                        return coord
                elif fmt == 'crystal':
                        return np.dot(self.vec,coord)*self.celldm
                elif fmt == 'alat':
                        return coord*self.celldm

        def get_coord(self,coord,fmt):
                if fmt == 'angstrom':
                        return coord*0.52917721092
                elif fmt == 'bohr':
                        return coord
                elif fmt == 'crystal':
                        return np.dot(self.vecinv,coord)/self.celldm
                elif fmt == 'alat':
                        return coord/self.celldm

        ######################################################
        # RETURN FUNCTIONS
        ######################################################
        def get_nat(self):
                return len(self.at_c)

        def get_celldm(self):
                return self.celldm

        def get_atom(self,index,fmt='bohr'):
                return [self.at_n[index],self.get_coord(self.at_c[index],fmt),fmt]

        def get_vec(self):
                return self.vec.T

        def get_comment(self):
                return self.comment

        def get_types(self):
                types = set()
                for i in self.at:
                    types.add(i.get_name())
                return types

        def get_ntyp(self):
                return len(self.get_types())

        def get_center(self):
                if not hasattr(self,'center'):
                        self.set_center()
                return self.center

        ######################################################
        # BOND FUNCTIONS
        ######################################################

        def get_bonds(self):
                if not hasattr(self,'bonds'): self.set_bonds()
                return self.bonds

        def set_bonds(self):
                nat = self.get_nat()
                self.bonds = []
                for i in range(nat):
                        for j in range(i,nat):
                                at_i = self.get_atom(i)
                                at_j = self.get_atom(j)
                                dist = np.linalg.norm(at_i[1]-at_j[1])
                                if at_i[0] != 'H' and at_j[0] != 'H':
                                        if 0.755 < dist < 3.5:
                                                self.bonds.append((i,j,dist))
                                else:
                                        if 0.755 < dist < 2.27:
                                                self.bonds.append((i,j,dist))



class PWParam(dict):

        def __init__(self):
                # make local copy of pse
                self['pse'] = copy.copy(pse)

                # k-point grid
                self['K_POINTS']=['gamma']
