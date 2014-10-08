#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import copy
from ptb_gui import MainWindow
from math import sqrt,floor
from time import clock
from collections import OrderedDict
import numpy as np
from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer

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
class TBController(QApplication):

        def __init__(self,argv):
                super(TBController,self).__init__(argv)
                self.argv = argv
                self.mol = []
                self.pwdata = []
                self.indict = OrderedDict([('xyz',self.parseXyz),
                               ('PWScf Input',self.parsePwi),
                               ('PWScf Output' , self.parsePwo),
                               ('PWO Final Conf.',self.parsePwoFinal)])
                self.outdict= OrderedDict([('PWScf Input',self.writePwi),
                               ('xyz',self.writeXyz)])
                QTimer.singleShot(0,self.argumentHandler)


#####################################################################
# Handle command line arguments:
#####################################################################
        def argumentHandler(self):
                #no argument: start GUI
                if len(self.argv) == 1:
                        self.gui = MainWindow(self)
                #check for misformatted options:
                elif self.argv[1][0]!='-': self.print_help()
                elif self.argv[1][0:2]=='--': self.print_help()
                #check for help request:
                elif '-h' in self.argv: self.print_help()
                #check for input files, start gui and load
                else:
                        self.gui = MainWindow(self)
                        #for i in range(1,len(self.argv)):
                        i=1
                        while i<len(self.argv):
                                if self.argv[i] == '-pwi':
                                        i+=1
                                        while i<len(self.argv) and self.argv[i][0]!='-':
                                                self.readFile('PWScf Input',self.argv[i])
                                                i+=1
                                        self.gui.centralWidget().loadView()
                                elif self.argv[i] == '-pwo':
                                        i+=1
                                        while i<len(self.argv) and self.argv[i][0]!='-':
                                                self.readFile('PWScf Output',self.argv[i])
                                                i+=1
                                        self.gui.centralWidget().loadView()
                                elif self.argv[i] == '-pwof':
                                        i+=1
                                        while i<len(self.argv) and self.argv[i][0]!='-':
                                                self.readFile('PWO Final Conf.',self.argv[i])
                                                i+=1
                                        self.gui.centralWidget().loadView()
                                elif self.argv[i] == '-xyz':
                                        i+=1
                                        while i<len(self.argv) and self.argv[i][0]!='-':
                                                self.readFile('xyz',self.argv[i])
                                                i+=1
                                        self.gui.centralWidget().loadView()

#####################################################################
# Print help
#####################################################################

        def print_help(self):
                f = sys.stdout
                f.write('PWToolBox usage:\n')
                f.write('ptb_main [OPTIONS]\n\n')
                f.write('No option given: start GUI\n\n')
                f.write('Options:\n')
                f.write('-h: print this help\n')
                f.write('-xyz [FILES]: open xyz file(s)\n')
                f.write('-pwi [FILES]: open PWScf input file(s)\n')
                f.write('-pwo [FILES]: open PWScf output file(s)\n')
                f.write('-pwof [FILES]: parse only the last config of PWO file(s)\n')
                self.quit()

#####################################################################
# Commandline Actions
#####################################################################

        #TODO: generate new files programmatically
        #def newFile(self,newfmt,newfile,oldcoord,oldparam=""):
        #def multFile(self,file,fmt,mult)

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
                                                tparam['K_POINTS']['active']='gamma'
                                        #MPGrid:
                                        #x y z offset
                                        #passed as whole string for now
                                        elif header[1] == 'automatic':
                                                line = data.pop(0).strip().split()
                                                tparam['K_POINTS']['automatic']=line
                                                tparam['K_POINTS']['active']='automatic'
                                        #else:
                                        #number of kpoints
                                        #x y z weight
                                        #passed as whole string for now
                                        else:
                                                nk = int(data.pop(0).strip().split()[0])
                                                kpoints = []
                                                for i in range(nk):
                                                        kpoints.append(data.pop(0).strip().split())
                                                        #kpoints.append(data.pop(0))
                                                #tparam['K_POINTS']=[header[1],nk,kpoints]
                                                tparam['K_POINTS']['disc']=kpoints
                                                tparam['K_POINTS']['active']=header[1]

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
                #TODO: temporarily saving everythin in ibrav=0 format, for gui-edit and saving
                #set celldm from input
                tmol.set_celldm(tparam['&system']['celldm(1)'])
                if tparam['&system']['ibrav'] == '0':
                        #check if CELL_PARAMETERS card has been read, if not present, throw error
                        if tvec == [[0,0,0],[0,0,0],[0,0,0]]:
                                print 'CELL_PARAMETERS missing'
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
                                print 'celldm(3) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                tmol.set_vec([[1,0,0],[-0.5,sqrt(3)*0.5,0],[0,0,ca]])
                elif tparam['&system']['ibrav'] == '5':
                        print 'cell descriptor not implemented yet'
                elif tparam['&system']['ibrav'] == '-5':
                        print 'cell descriptor not implemented yet'
                elif tparam['&system']['ibrav'] == '6':
                        #simple tetragonal
                        if not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                tmol.set_vec([[1,0,0],[0,1,0],[0,0,ca]])
                elif tparam['&system']['ibrav'] == '7':
                        #body centered tetragonal
                        if not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                tmol.set_vec([[0.5,-0.5,ca*0.5],[0.5,0.5,ca*0.5],[-0.5,-0.5,ca*0.5]])
                elif tparam['&system']['ibrav'] == '8':
                        #simple orthorhombic
                        if not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                tmol.set_vec([[1,0,0],[0,ba,0],[0,0,ca]])
                elif tparam['&system']['ibrav'] == '9':
                        #basis centered orthorhombic
                        if not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                tmol.set_vec([[0.5,ba*0.5,0],[-0.5,ba*0.5,0],[0,0,ca]])
                elif tparam['&system']['ibrav'] == '10':
                        #face centered orthorhombic
                        if not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                tmol.set_vec([[0.5,0,ca*0.5],[0.5,ba*0.5,0],[0,ba*0.5,ca*0.5]])
                elif tparam['&system']['ibrav'] == '11':
                        #body centered orthorhombic
                        if not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                tmol.set_vec([[0.5,ba*0.5,ca*0.5],[-0.5,ba*0.5,ca*0.5],[-0.5,-ba*0.5,ca*0.5]])
                elif tparam['&system']['ibrav'] == '12':
                        #simple monoclinic
                        if not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        elif not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(4)' in tparam['&system']:
                                print 'celldm(4) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                cg = float(tparam['&system']['celldm(4)'])
                                tmol.set_vec([[1,0,0],[ba*cg,ba*sqrt(1-cg),0],[0,0,ca]])
                elif tparam['&system']['ibrav'] == '-12':
                        #simple monoclinic, alternate definition
                        if not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        elif not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(5)' in tparam['&system']:
                                print 'celldm(5) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                cb = float(tparam['&system']['celldm(5)'])
                                tmol.set_vec([[1,0,0],[0,ba,0],[ca*cb,0,ca*sqrt(1-cb)]])
                elif tparam['&system']['ibrav'] == '13':
                        #base centered monoclinic
                        if not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        elif not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(4)' in tparam['&system']:
                                print 'celldm(4) missing'
                        else:
                                ca = float(tparam['&system']['celldm(3)'])
                                ba = float(tparam['&system']['celldm(2)'])
                                cg = float(tparam['&system']['celldm(4)'])
                                tmol.set_vec([[0.5,0,-ca*0.5],[ba*cg,ba*sqrt(1-cg),0],[0.5,0,ca*0.5]])
                elif tparam['&system']['ibrav'] == '14':
                        #base centered monoclinic
                        if not 'celldm(2)' in tparam['&system']:
                                print 'celldm(2) missing'
                        elif not 'celldm(3)' in tparam['&system']:
                                print 'celldm(3) missing'
                        elif not 'celldm(4)' in tparam['&system']:
                                print 'celldm(4) missing'
                        elif not 'celldm(5)' in tparam['&system']:
                                print 'celldm(5) missing'
                        elif not 'celldm(6)' in tparam['&system']:
                                print 'celldm(6) missing'
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
                #TODO: recreate used parameters
                #tparam=PWParam()
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
                                tmol = Molecule()
                                tmol.set_celldm(celldm)
                                tmol.set_vec(vec)
                                for j in range(i+1,i+nat+1):
                                        atom = data[j].split()
                                        tmol.create_atom(atom[1],float(atom[6]),float(atom[7]),float(atom[8]),'alat')
                                i+=nat
                                tlist.append(tmol)
                        #read step-vectors if cell is variable
                        elif line[0] == 'CELL_PARAMETERS':
                                for j in [0,1,2]:
                                        temp = data[i+1+j].split()
                                        vec[j]=[float(x) for x in temp[0:3]]
                        #read step-coordinates
                        elif line[0] == 'ATOMIC_POSITIONS':
                                tmol = Molecule()
                                tmol.set_celldm(celldm)
                                tmol.set_vec(vec)
                                for j in range(i+1,i+nat+1):
                                        atom = data[j].split()
                                        tmol.create_atom(atom[0],float(atom[1]),float(atom[2]),float(atom[3]),line[1].strip('()'))
                                i+=nat
                                tlist.append(tmol)
                        #break on reaching final coordinates (duplicate)
                        elif line[0] == 'Begin':
                                break
                        #ignore everything else
                        else:
                                pass
                        i+=1
                self.mol.append(tlist)

        def parsePwoFinal(self,data):
                #parse only the final config for commandline actions
                i=0
                vec=[[0,0,0],[0,0,0],[0,0,0]]
                final = False
                while i<len(data):
                        line = data[i].split()
                        #print line
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
                                tmol = Molecule()
                                tmol.set_celldm(celldm)
                                tmol.set_vec(vec)
                                for j in range(i+1,i+nat+1):
                                        atom = data[j].split()
                                        tmol.create_atom(atom[0],float(atom[1]),float(atom[2]),float(atom[3]),line[1].strip('()'))
                        else:
                                pass
                        i+=1
                self.mol.append([tmol])

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
                        atom=mol.get_atom(i,coordfmt)
                        f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                            atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
                f.write('\n')

                #K_POINTS
                f.write('K_POINTS'+' '+param['K_POINTS']['active']+'\n')
                #Gamma point only
                if param['K_POINTS']['active'] == 'gamma':
                        pass
                #MPGrid:
                #x y z offset
                elif param['K_POINTS']['active'] == 'automatic':
                        f.write('{:4s}{:4s}{:4s}{:4d}{:4d}{:4d}'.format(
                                param['K_POINTS']['automatic'][0],
                                param['K_POINTS']['automatic'][1],
                                param['K_POINTS']['automatic'][2],
                                param['K_POINTS']['automatic'][3],
                                param['K_POINTS']['automatic'][4],
                                param['K_POINTS']['automatic'][5]
                                )+'\n')
                #number of kpoints
                #x y z weight
                else:
                        f.write(str(len(param['K_POINTS']['disc']))+'\n')
                        for i in range(len(param['K_POINTS']['disc'])):
                                f.write('{:4s}{:4s}{:4s}{:4s}'.format(
                                        param['K_POINTS']['disc'][i][0],
                                        param['K_POINTS']['disc'][i][1],
                                        param['K_POINTS']['disc'][i][2],
                                        param['K_POINTS']['disc'][i][3]
                                        )+'\n')
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
                self.at_n.append(name)
                self.at_c.append(self.set_coord([x,y,z],fmt))

        # append copy of existing atom
        def append_atom_cp(self,addat):
                self.at_n.append(addat[0])
                self.at_c.append(self.set_coord(addat[1],addat[2]))

        # insert atom at given position
        def insert_atom(self,pos,addat):
                self.at_n.insert(pos,addat[0])
                self.at_c.insert(pos,self.set_coord(addat[1],addat[2]))

        # remove atom
        def del_atom(self,index):
                del self.at_n[index]
                del self.at_c[index]

        # append molecule
        #def append_mol(self, mol):
        #        for i in range(mol.get_nat()):
        #                self.at.append(copy.copy(mol.at[i]))

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

        # set vectors in bohr
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
                for i in self.at_n:
                    types.add(i)
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

        def get_pbc_bonds(self):
                if not hasattr(self,'pbc_bonds'): self.set_pbc_bonds()
                return self.pbc_bonds

        def set_bonds(self):
                nat = self.get_nat()
                self.bonds = []
                at_c = self.at_c
                at_n = self.at_n

                for i in range(nat):
                        for j in range(i+1,nat):
                                dist = at_c[i]-at_c[j]

                                #cancel if distance in one direction is greater than allowed bond length
                                if dist[0]>3.5 or dist[1]>3.5 or dist[2]>3.5: continue

                                dist = np.dot(dist,dist)
                                if at_n[i] != 'H' and at_n[j] != 'H':
                                        #maximum bond length: 1.9A
                                        if 0.57 < dist < 12.25:
                                                self.bonds.append([at_c[i],at_c[j],at_n[i],at_n[j]])
                                else:
                                        #maximum bond length for hydrogen: 1.2A
                                        if 0.57 < dist < 5.15:
                                                self.bonds.append([at_c[i],at_c[j],at_n[i],at_n[j]])

        def set_pbc_bonds(self):
                nat = self.get_nat()
                self.pbc_bonds=[self.get_bonds(),[],[],[],[],[],[],[]]
                vec = self.get_vec()*self.get_celldm()
                off = [0,vec[0],vec[1],vec[2],vec[0]+vec[1],vec[0]+vec[2],vec[1]+vec[2],vec[0]+vec[1]+vec[2]]
                at_c = self.at_c
                at_n = self.at_n
                for i in range(nat):
                        for j in range(nat):
                                dist_at = self.at_c[i] - self.at_c[j]
                                for k in [1,2,3,4,5,6,7]:
                                        dist = dist_at+off[k]

                                        #cancel if distance in one direction is greater than allowed bond length
                                        if dist[0]>3.5 or dist[1]>3.5 or dist[2]>3.5: continue

                                        dist = np.dot(dist,dist)
                                        if at_n[i] != 'H' and at_n[j] != 'H':
                                                #maximum bond length: 1.9A
                                                if 0.57 < dist < 12.25:
                                                        self.pbc_bonds[k].append([at_c[i]+off[k],at_c[j],at_n[i],at_n[i]])
                                        else:
                                                #maximum bond length for hydrogen: 1.2A
                                                if 0.57 < dist < 5.15:
                                                        self.pbc_bonds[k].append([at_c[i]+off[k],at_c[j],at_n[i],at_n[i]])

        #####################################################
        # EDIT FUNCTIONS
        #####################################################

        def mult(self,x,y,z):
                nat = self.get_nat()
                vec = self.get_vec()*self.get_celldm()
                mult = [x,y,z]
                for k in [0,1,2]:
                        for i in range(1,mult[k]):
                                for j in range(nat):
                                        self.append_atom_cp(self.get_atom(j))
                                        atom = self.get_atom(-1)
                                        self.set_atom(-1,atom[0],atom[1]+i*vec[k],'bohr')
                        nat = self.get_nat()
                self.set_vec(self.get_vec()*[[x],[y],[z]])


class PWParam(dict):

        def __init__(self):
                # make local copy of pse
                self['pse'] = copy.copy(pse)

                # k-point grid
                self['K_POINTS']={'active':'gamma'}
