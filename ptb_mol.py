#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import copy

######################################################################
# PSE DICTIONARY ? suitable for lammps ?
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

        def readFile(self,fmt,filename):
                parser = {'PWi' : self.parsePwi,
                          'PWo' : self.parsePwo,
                          'xyz' : self.parseXyz,
                         }
                data = open(filename,'r').readlines()
                parser[fmt](data)
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
                        nat = int(data[0].strip())
                        tmol.comment = data[1].strip()
                        # read coordinates and types
                        for i in range(2,nat+2):
                                line = data[i].split()
                                tmol.create_atom(line[0],float(line[1]),float(line[2]),float(line[3]))
                        # remove parsed lines to check if there's more to do
                        for i in range(0,nat+2):
                                del data[0]
                        tlist.append(copy.deepcopy(tmol))
                        del tmol
                #return tlist
                self.mol = self.mol + tlist

        def parsePwi(self,data):
                # no need for list, only one molecule per file
                tmol = Molecule()
                tparam = PWParam()
                k = 0
                while data:
                        k +=1
                        # ignore empty lines
                        if not data[0].strip():
                                del data[0]
                                print 'empty line'
                                print data
                        # parse namelists
                        elif data[0].strip()[0] == '&':
                                tnl = tparam.Namelist()
                                i = 1
                                # parse entries
                                while data[i].strip() != '/':
                                        line = data[i].strip().split(',')
                                        for j in range(len(line)):
                                                tnl[line[j].split('=')[0]]=line[j].split('=')[1]
                                        i +=1
                                tparam[data[0].strip()]=tnl
                                # delete parsed namelist
                                for j in range(i+1):
                                        del data[0]
                                print 'new namelist'
                                print data
                        # parse card
                        elif data[0].strip().split()[0].isupper():
                                i = 1
                                # 7 types of cards, need hardcoding
                                # 4 supported for now
                                # still need modes for all of them
                                while data[i].strip():
                                        line = data[i].strip().split()
                                        if data[0].strip().split()[0] == 'ATOMIC_SPECIES':
                                                tparam.pse[line[0]][1] = float(line[1])
                                                tparam.pse[line[0]].append(line[2])
                                        elif data[0].strip().split()[0] == 'ATOMIC_POSITIONS':
                                                tmol.create_atom(line[0],float(line[1]),float(line[2]),float(line[3]))
                                        elif data[0].strip().split()[0] == 'K_POINTS':
                                                tparam['K_POINTS']=[line[0:3],line[3:7]]
                                        elif data[0].strip().split()[0] == 'CELL_PARAMETERS':
                                                tmol.vec[i-1]=[float(line[0]),float(line[1]),float(line[2])]
                                        i +=1
                                # delete parsed card
                                for j in range(i):
                                        del data[0]
                                print 'new card'
                                print data
                print k
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
                self.vec=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
                self.offset=[0.0,0.0,0.0]
                self.comment = ''

        # append new atom
        def create_atom(self,name,x,y,z):
                self.at.append(self.Atom(name,x,y,z))
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

        # set periodicity
        def set_periodicity(self,vec0,vec1,vec2,off=[0.0,0.0,0.0]):
                self.vec[0]=vec0
                self.vec[1]=vec1
                self.vec[2]=vec2
                self.offset=off

######################################################################
# ATOM CLASS
######################################################################

        class Atom:

                def __init__(self,name,x=0,y=0,z=0):
                        self.name=name
                        self.number=pse[self.name][0]
                        self.weight=pse[self.name][1]
                        self.coord=[x,y,z]

                ######################################################
                # return functions
                ######################################################
                def get_name(self):
                        return self.name

                def get_number(self):
                        return self.number

                def get_weight(self):
                        return self.weight

                def get_coord(self):
                        return self.coord

                #############################################################
                # set functions
                #############################################################
                def set(self,name,x,y,z):
                        self.name=name
                        self.number=pse[name][0]
                        self.weight=pse[name][1]
                        self.coord=[x,y,z]

        #Kommt spaeter:
        #class Bond:

class PWParam(dict):

        def __init__(self):
                # make local copy of pse
                self.pse = copy.copy(pse)

                # k-point grid
                # only MP-grids (automatic) for now
                self['K_POINTS']=[[1,1,1],[0,0,0]]

        # define NL and Card classes just to distinguish
        class Namelist(dict):
                def is_namelist(self):
                        return True

        #class Card(dict):
        #        def is_namelist(self):
        #                return False
