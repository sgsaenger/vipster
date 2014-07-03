#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import copy
import pwtool
from math import sqrt
import numpy as np
from PyQt4.QtGui import *

######################################################################
# PSE DICTIONARY
######################################################################
# pse[0] id
# pse[1] ~ weight
pse={"X":  [0, 0.0],

     "H" : [1, 1.0079],
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

        def __init__(self):
                self.mol = []
                self.pwdata = []
                #self.data = []
                self.gui = pwtool.CoordTB(self)
                self.indict = {'xyz' : self.parseXyz,
                               'PWScf Input' : self.parsePwi,
                               'PWScf Output' : self.parsePwo,
                               }
                self.outdict= {'PWScf Input' : self.parsePwi,
                               'xyz' : self.parseXyz,
                               }

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
                while data:
                        # handle empty lines at eof or between molecules
                        if not data[0].strip().isdigit():
                                del data[0]
                                if not data: break
                        # create new molecule
                        tmol = Molecule()
                        # fixed format nat and comment
                        nat = int(data.pop(0).strip())
                        tmol.comment = data.pop(0).strip()
                        # read coordinates and types
                        for i in range(nat):
                                line = data.pop(0).split()
                                tmol.create_atom(line[0],float(line[1]),float(line[2]),float(line[3]))
                        # remove parsed lines to check if there's more to do
                        #for i in range(0,nat+2):
                        #        del data[0]
                        tlist.append(tmol)
                #return tlist
                self.mol = self.mol + tlist

        def parsePwi(self,data):
                # no need for list, only one molecule per file
                tmol = Molecule()
                tparam = PWParam()
                while data:
                        header = data.pop(0).strip().split()
                        # ignore empty lines
                        if not header:
                                #debug output. case not really needed.
                                print 'empty line'
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
                                print 'new namelist', header[0]
                                print tparam
                        # parse card
                        elif header[0][0].isupper():
                                # 7 types of cards, need hardcoding
                                # 4 supported for now
                                # still need modes for all of them
                                if header[0] == 'ATOMIC_SPECIES':
                                        for i in range(int(tparam['&system']['ntyp'])):
                                                line = data.pop(0).strip().split()
                                                tparam['pse'][line[0]][1] = float(line[1])
                                                tparam['pse'][line[0]].append(line[2])
                                elif header[0] == 'ATOMIC_POSITIONS':
                                        for i in range(int(tparam['&system']['nat'])):
                                                line = data.pop(0).strip().split()
                                                tmol.create_atom(line[0],float(line[1]),float(line[2]),float(line[3]))
                                elif header[0] == 'K_POINTS':
                                        if header[1] == 'automatic':
                                                line = data.pop(0).strip().split()
                                                tparam['K_POINTS']=[line[0:3],line[3:7]]
                                elif header[0] == 'CELL_PARAMETERS':
                                        vec=[[0,0,0],[0,0,0],[0,0,0]]
                                        for i in range(3):
                                                line = data.pop(0).strip().split()
                                                vec[i]=[float(x)for x in line]
                                        tmol.set_vec(vec)
                                        print tmol.get_vec()
                                else:
                                        print 'unknown or unsupported card.'
                                print 'new card', header[0]
                self.mol.append(tmol)
                self.pwdata.append(tparam)

        def parsePwo(self,data):
                bargl = 0

        def writeXyz(self,mol,filename=""):
                if filename == "":
                        f=sys.stdout
                else:
                        f=open(filename,'w')
                # fixed format nat and comment
                f.write(str(mol.get_nat())+'\n')
                f.write(mol.get_comment()+'\n')
                # write coordinates
                for cntat in range(0,mol.get_nat()):
                        f.write(
                                '{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                                     mol.at[cntat].get_name(),
                                     mol.at[cntat].get_coord()[0],
                                     mol.at[cntat].get_coord()[1],
                                     mol.at[cntat].get_coord()[2]
                                     )+'\n'
                                )

        def writePwi(self,mol,param,coordfmt='Crystal',filename=""):
                if filename == "":
                        f=sys.stdout
                else:
                        f=open(filename,'w')
                #
                

######################################################################
# MOLECULE CLASS
######################################################################
class Molecule:

        def __init__(self):
                # set molecule info
                #self.nat = 0
                # set atom list
                self.at=[]
                self.celldm = 1.0
                self.vec=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self.vecinv=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self.offset=[0.0,0.0,0.0]
                self.comment = ''

        # append new atom
        def create_atom(self,name,x,y,z):
                self.at.append(self.Atom(self,name,x,y,z))
                #self.nat = self.nat + 1

        # append copy of existing atom
        def append_atom_cp(self,addat):
                self.at.append(copy.copy(addat))
                #self.nat = self.nat + 1

        # append molecule
        def append_mol(self, mol):
                for i in range(mol.get_nat()):
                        self.at.append(copy.copy(mol.at[i]))
                #self.nat = self.nat + mol.nat

        # remove atom
        def del_atom(self,index):
                del self.at[index]

        ######################################################
        # return functions
        ######################################################
        def get_nat(self):
                #return self.nat
                return len(self.at)

        def get_celldm(self):
                return self.celldm

        def get_atom(self,atom):
                return self.at[atom]

        def get_vec(self):
                return self.vec.T

        def get_comment(self):
                return self.comment

        def get_ntype(self):
                types = set()
                for i in self.at:
                    types.add(i.get_name())
                return len(types)

        ######################################################
        # set functions
        ######################################################

        # set celldm
        def set_celldm(self,cdm):
                self.celldm = cdm

        # set vectors
        def set_vec(self,vec):
                self.vec = np.array(vec).T

        def set_periodicity(self,vec0,vec1,vec2,off=[0.0,0.0,0.0]):
                self.vec = np.array([vec0,vec1,vec2]).T
                self.vecinv = np.linalg.inv(self.vec)
                self.offset=off

######################################################################
# ATOM CLASS
######################################################################

        class Atom:

                def __init__(self,mol,name,x=0,y=0,z=0):
                        self.name=name
                        self.number=pse[self.name][0]
                        self.weight=pse[self.name][1]
                        self.coord=np.array([x,y,z])
                        self.mol = mol

                ######################################################
                # return functions
                ######################################################
                def get_name(self):
                        return self.name

                def get_number(self):
                        return self.number

                def get_weight(self):
                        return self.weight

                def get_coord(self,fmt):
                        if fmt == u'Ångström':
                                return self.coord
                        elif fmt == 'Bohr':
                                return self.coord/0.52917721092
                        elif fmt == 'Crystal':
                                return np.dot(self.mol.vecinv,self.coord)/self.mol.celldm
                        elif fmt == 'Alat':
                                return self.coord/self.mol.celldm

                #############################################################
                # set functions
                #############################################################
                def set(self,name,x,y,z):
                        self.name=name
                        self.number=pse[name][0]
                        self.weight=pse[name][1]
                        self.coord=np.array([x,y,z])

                def set_name(self,name):
                        self.name = name

                def set_coord(self,fmt,coord):
                        self.coord=np.array(coord)
                        if fmt == u'Ångström':
                                pass
                        elif fmt == 'Bohr':
                                self.coord *= 0.52917721092
                        elif fmt == 'Crystal':
                                self.coord = np.dot(self.mol.vec,self.coord)*self.mol.celldm
                        elif fmt == 'Alat':
                                self.coord *= self.mol.celldm


        #Kommt spaeter:
        #class Bond:

class PWParam(dict):

        def __init__(self):
                # make local copy of pse
                self['pse'] = copy.copy(pse)

                # k-point grid
                # only MP-grids (automatic) for now
                self['K_POINTS']=[[1,1,1],[0,0,0]]

        # define NL and Card classes just to distinguish
        #class Namelist(dict):
        #        def is_namelist(self):
        #                return True

        #class Card(dict):
        #        def is_namelist(self):
        #                return False

#####################################################
# Application
#####################################################

def main():
        app = QApplication(sys.argv)
        control = TBController()
        sys.exit(app.exec_())

if __name__ == '__main__':
        main()
