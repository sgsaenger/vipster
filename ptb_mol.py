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
                                                tcoord.append(data.pop(0).strip().split())

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

                #Append to controller
                self.mol.append(tmol)
                self.pwdata.append(tparam)

        def parsePwo(self,data):
                #Multiple configs supported
                tlist = []
                #no parameters for now
                #tparam=PWParam()

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
                        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                                     mol.at[cntat].get_name(),
                                     mol.at[cntat].get_coord('Bohr')[0],
                                     mol.at[cntat].get_coord('Bohr')[1],
                                     mol.at[cntat].get_coord('Bohr')[2]
                                     )+'\n'
                                )
                f.close()

        def writePwi(self,mol,param,coordfmt='crystal',filename=""):
                if filename == "":
                        f=sys.stdout
                else:
                        f=open(filename,'w')

                #set nat and ntyp

                #&control, &system and &electron namelists are mandatory
                for i in ['&control','&system','&electrons']:
                        f.write(i+'\n')
                        #write all parameters
                        for j in range(len(param[i])):
                            f.write(param[i].keys()[j]+'='+param[i].values()[j]+'\n')
                        f.write('/\n\n')
                #&ions only when needed
                if param['&control']['calculation'] in ['relax','vc-relax','md','vc-md']:
                        f.write('&ions'+'\n')
                        for j in range(len(param['&ions'])):
                                f.write(param['&ions'].keys()[j]+'='+param['&ions'].values()[j]+'\n')
                        f.write('/\n\n')
                #&cell only when needed
                if param['&control']['calculation'] in ['vc-relax','vc-md']:
                        f.write('&cell'+'\n')
                        for j in range(len(param['&cell'])):
                                f.write(param['&cell'].keys()[j]+'='+param['&cell'].values()[j]+'\n')
                        f.write('/\n\n')

                #ATOMIC_SPECIES card:
                f.write('ATOMIC_SPECIES'+'\n')
                types = list(mol.get_types())
                for i in range(len(mol.get_types())):
                        atom = types[i]
                        f.write(atom+'    '+str(param['pse'][atom][1])+'   '+str(param['pse'][atom][2])+'\n')
                f.write('\n')

                #ATOMIC_POSITIONS
                f.write('ATOMIC_SPECIES'+' '+coordfmt+'\n')
                for i in range(mol.get_nat()):
                        atom=mol.get_atom(i)
                        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                            atom.get_name(),
                            atom.get_coord(coordfmt)[0],
                            atom.get_coord(coordfmt)[1],
                            atom.get_coord(coordfmt)[2],
                            )+'\n'
                            )
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
                self.at=[]
                self.celldm = 1.0
                self.vec=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self.vecinv=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self.offset=[0.0,0.0,0.0]
                self.comment = ''

        # append new atom
        def create_atom(self,name,x=0.,y=0.,z=0.,fmt='angstrom'):
                self.at.append(self.Atom(self,name,x,y,z,fmt))
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
                return len(self.at)

        def get_celldm(self):
                return self.celldm

        def get_atom(self,atom):
                return self.at[atom]

        def get_vec(self):
                return self.vec.T

        def get_comment(self):
                return self.comment

        def get_types(self):
                types = set()
                for i in self.at:
                    types.add(i.get_name())
                return types

        ######################################################
        # set functions
        ######################################################

        # set celldm
        def set_celldm(self,cdm):
                self.celldm = float(cdm)

        # set vectors
        def set_vec(self,vec):
                self.vec = np.array(vec).T
                self.vecinv = np.linalg.inv(self.vec)

        #def set_periodicity(self,vec0,vec1,vec2,off=[0.0,0.0,0.0]):
        #        self.vec = np.array([vec0,vec1,vec2]).T
        #        self.vecinv = np.linalg.inv(self.vec)
        #        self.offset=off

######################################################################
# ATOM CLASS
######################################################################

        class Atom:

                def __init__(self,mol,name,x,y,z,fmt):
                        self.name=name
                        #self.number=pse[self.name][0]
                        #self.weight=pse[self.name][1]
                        self.mol = mol
                        self.set_coord(fmt,[x,y,z])

                ######################################################
                # return functions
                ######################################################
                def get_name(self):
                        return self.name

                #def get_number(self):
                #        return self.number

                #def get_weight(self):
                #        return self.weight

                def get_coord(self,fmt):
                        if fmt == 'angstrom':
                                return self.coord*0.52917721092
                        elif fmt == 'bohr':
                                return self.coord
                        elif fmt == 'crystal':
                                return np.dot(self.mol.vecinv,self.coord)/self.mol.celldm
                        elif fmt == 'alat':
                                return self.coord/self.mol.celldm

                #############################################################
                # set functions
                #############################################################
                #def set(self,name,x,y,z):
                #        self.name=name
                #        self.number=pse[name][0]
                #        self.weight=pse[name][1]
                #        self.coord=np.array([x,y,z])

                def set_name(self,name):
                        self.name = name

                def set_coord(self,fmt,coord):
                        self.coord=np.array(coord)
                        if fmt == 'angstrom':
                                self.coord /=0.52917721092
                        elif fmt == 'bohr':
                                #self.coord *= 0.52917721092
                                pass
                        elif fmt == 'crystal':
                                self.coord = np.dot(self.mol.vec,self.coord)*self.mol.celldm
                        elif fmt == 'alat':
                                self.coord *= self.mol.celldm


        #Kommt spaeter:
        #class Bond:

class PWParam(dict):

        def __init__(self):
                # make local copy of pse
                self['pse'] = copy.copy(pse)

                # k-point grid
                self['K_POINTS']=['gamma']

                # create dummy namelists
                self['&control']=[]
                self['&system']=[]
                self['&electrons']=[]
                self['&ions']=[]
                self['&cell']=[]

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
