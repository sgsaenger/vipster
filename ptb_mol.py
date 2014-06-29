#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import copy

######################################################################
# MAIN CONTROLLER CLASS
######################################################################
class TBController:

        def __init__(self):
                self.mol = []
                self.pwdata = []
                #self.data = []

        def readFile(self,fmt,filename):
                parser = {'PWi' : self.readPwi,
                          'PWo' : self.readPwo,
                          'xyz' : self.parseXyz,
                         }
                data = open(filename,'r').readlines()
                parser[fmt](data)
                #return data

        def parseXyz(self,data):
                tlist = []
                molcnt = 0
                if not data[0].strip().isdigit():
                        print 'not a xyz file'
                        return
                while data:
                        molcnt += 1
                        if not data[0].strip().isdigit():
                                del data[0]
                                if not data: break
                        tmol = Molecule()
                        nat = int(data[0].strip())
                        tmol.comment = data[1].strip()
                        for i in range(2,nat+2):
                                line = data[i].split()
                                tmol.create_atom(line[0],float(line[1]),float(line[2]),float(line[3]))
                        for i in range(0,nat+2):
                                del data[0]
                        tlist.append(copy.deepcopy(tmol))
                        del tmol
                #return tlist
                self.mol = self.mol + tlist

        def readPwi(self):
                bargl = 0

        def readPwo(self):
                bargl = 0

        def writeXyz(self,mol,filename=""):
                if filename == "":
                        f=sys.stdout
                else:
                        f=open(filename,'w')
                f.write(str(mol.get_nat())+'\n')
                f.write(mol.get_comment()+'\n')
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

        # set periodicity
        def set_periodicity(self,vec0,vec1,vec2,off=[0.0,0.0,0.0]):
                self.vec[0]=vec0
                self.vec[1]=vec1
                self.vec[2]=vec2
                self.offset=off

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

class PWParam:

        def __init__(self):
                self.calctype = 'scf'
                self.ntyp = 0


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
