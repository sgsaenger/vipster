#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import copy
from math import sqrt
from collections import OrderedDict
from PyQt4.QtGui import *
from PyQt4.QtCore import QTimer

from molecule import Molecule
from ptb_gui import MainWindow

######################################################################
# PSE DICTIONARY
######################################################################
# pse[0] id
# pse[1] ~ weight
pse={"X":  [0,0.0],

     "H" : [1,1.0079,'H.uspp736.pbe.UPF'],
     "He": [2,4.0026,'He.uspp736.pbe.UPF'],

     "Li": [3, 6.941],
     "Be": [4, 9.0122],
     "B" : [5,10.811,'B.uspp736.pbe.UPF'],
     "C" : [6,12.0107,'C.uspp736.pbe.UPF'],
     "N" : [7,14.007,'N.uspp736.pbe.UPF'],
     "O" : [8,15.999,'O.uspp736.pbe.UPF'],
     "F" : [9,18.998,'F.uspp736.pbe.UPF'],
     "Ne": [10,20.18],

     "Na": [11,22.99],
     "Mg": [12,24.305],
     "Al": [13,26.982,'Al.uspp736.pbe.UPF'],
     "Si": [14,28.086,'Si.uspp736.pbe.UPF'],
     "P" : [15,30.974,'P.uspp736.pbe.UPF'],
     "S" : [16,32.065,'S.uspp736.pbe.UPF'],
     "Cl": [17,35.453,'Cl.uspp736.pbe.UPF'],
     "Ar": [18,39.948],

     "K" : [19,39.098],
     "Ca": [20,40.078],
     "Sc": [21,44.9559],
     "Ti": [22,47.867],
     "V" : [23,50.9415],
     "Cr": [24,51.9961],
     "Mn": [25,54.938],
     "Fe": [26,55.845],
     "Co": [27,58.9332],
     "Ni": [28,58.6934],
     "Cu": [29,63.546],
     "Zn": [30,65.39],
     "Ga": [31,69.723],
     "Ge": [32,72.64],
     "As": [33,74.922],
     "Se": [34,78.96],
     "Br": [35,79.904,'Br.uspp736.pbe.UPF'],
     "Kr": [36,83.798],

     "Rb" : [37,85.4678],
     "Sr" : [38,87.62],
     "Y"  : [39,88.9059],
     "Zr" : [40,91.224],
     "Nb" : [41,92.9064],
     "Mo" : [42,95.94],
     "Tc" : [43,98.0],
     "Ru" : [44,101.07],
     "Rh" : [45,102.9055],
     "Pd" : [46,106.42],
     "Ag" : [47,107.8682],
     "Cd" : [48,112.411],
     "In" : [49,114.818],
     "Sn" : [50,118.71],
     "Sb" : [51,121.76],
     "Te" : [52,127.6],
     "I"  : [53,126.9045],
     "Xe" : [54,131.293],

     "Cs" : [55,132.9055],
     "Ba" : [56,137.327],
     "La" : [57,138.9055],
     "Ce" : [58,140.116],
     "Pr" : [59,140.9077],
     "Nd" : [60,144.24],
     "Pm" : [61,145.0],
     "Sm" : [62,150.36],
     "Eu" : [63,151.964],
     "Gd" : [64,157.25],
     "Tb" : [65,158.9253],
     "Dy" : [66,162.5],
     "Ho" : [67,164.9303],
     "Er" : [68,167.259],
     "Tm" : [69,168.9342],
     "Yb" : [70,173.04],
     "Lu" : [71,174.967],
     "Hf" : [72,178.49],
     "Ta" : [73,180.9479],
     "W"  : [74,183.84],
     "Re" : [75,186.207],
     "Os" : [76,190.23],
     "Ir" : [77,192.217],
     "Pt" : [78,195.078],
     "Au" : [79,196.9665],
     "Hg" : [80,200.592],
     "Tl" : [81,204.3833],
     "Pb" : [82,207.2],
     "Bi" : [83,208.9804],
     "Po" : [84,209],
     "At" : [85,210],
     "Rn" : [86,222],

     "Fr" : [87,223],
     "Ra" : [88,226],
     "Ac" : [89,227],
     "Th" : [90,232.03806],
     "Pa" : [91,231.03588],
     "U"  : [92,238.02891],
     "Np" : [93,237],
     "Pu" : [94,244],
     "Am" : [95,243],
     "Cm" : [96,247],
     "Bk" : [97,247],
     "Cf" : [98,251],
     "Es" : [99,252],
     "Fm" : [100,257],
     "Md" : [101,258],
     "No" : [102,259],
     "Lr" : [103,266],
     "Rf" : [104,267],
     "Db" : [105,268],
     "Sg" : [106,269],
     "Bh" : [107,270],
     "Hs" : [108,269],
     "Mt" : [109,278],
     "Ds" : [110,281],
     "Rg" : [111,281],
     "Cn" : [112,285],
     "Uut": [113,286],
     "Fl" : [114,289],
     "Uup": [115,289],
     "Lv" : [116,293],
     "Uus": [117,294],
     "Uuo": [118,294]}

######################################################################
# MAIN CONTROLLER CLASS
######################################################################
class TBController(QApplication):

        def __init__(self,argv):
                super(TBController,self).__init__(argv)
                self.argv = argv
                self._mol = []
                self._pwdata = []
                self.indict = OrderedDict([('xyz',self._parseXyz),
                               ('PWScf Input',self._parsePwi),
                               ('PWScf Output' , self._parsePwo),
                               ('PWO Final Conf.',self._parsePwoFinal),
                               ('Gaussian Cube File',self._parseCube)])
                self.outdict= OrderedDict([('PWScf Input',self._writePwi),
                               ('xyz',self._writeXyz)])
                QTimer.singleShot(0,self._argumentHandler)


#####################################################################
# Handle command line arguments:
#####################################################################
        def _argumentHandler(self):
                #no argument: start GUI
                if len(self.argv) == 1:
                        self.gui = MainWindow(self)
                #check for misformatted options or help request:
                elif self.argv[1][0]!='-': self._print_help(22)
                elif self.argv[1][0:2]=='--': self._print_help(22)
                elif '-h' in self.argv: self._print_help(0)
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
                                elif self.argv[i] == '-cube':
                                        i+=1
                                        while i<len(self.argv) and self.argv[i][0]!='-':
                                                self.readFile('Gaussian Cube File',self.argv[i])
                                                i+=1
                                        self.gui.centralWidget().loadView()

#####################################################################
# Print help
#####################################################################

        def _print_help(self,err):
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
                f.write('-cube [FILES]: open CUBE file(s)\n')
                self.exit(err)

#####################################################################
# GET FUNCTIONS
#####################################################################

        def get_mol(self,index,step):
                return self._mol[index][step]

        def get_lmol(self,index):
                return len(self._mol[index])

        def get_nmol(self):
                return len(self._mol)

        def get_pw(self,index):
                return self._pwdata[index]

        def get_npw(self):
                return len(self._pwdata)

#####################################################################
# NEW MOLECULE
#####################################################################

	def newMol(self):
		self._mol.append([Molecule()])
#####################################################################
# READ FUNCTIONS
#####################################################################

        def readFile(self,fmt,filename):
                data = open(filename,'r').readlines()
                self.indict[fmt](data)
                #return data

        def _parseXyz(self,data):
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
                self._mol.append(tlist)

        def _parsePwi(self,data):
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
                        #trigonal
                        if not 'celldm(4)' in tparam['&system']:
                                print 'celldm(4) missing'
                        else:
                                c = float(tparam['&system']['celldm(4)'])
                                tx=sqrt((1-c)/2)
                                ty=sqrt((1-c)/6)
                                tz=sqrt((1+2*c)/3)
                                tmol.set_vec([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]])
                elif tparam['&system']['ibrav'] == '-5':
                        #trigonal,alternative
                        if not 'celldm(4)' in tparam['&system']:
                                print 'celldm(4) missing'
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
                self._mol.append(tlist)

        def _parsePwoFinal(self,data):
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
                self._mol.append([tmol])

        def _parseCube(self,data):
                tmol = Molecule()
                tcoord=[]
                tvec=[[0,0,0],[0,0,0],[0,0,0]]
                #parse data
                i=0
                #two lines of comments, combine
                tmol.comment=data[0]+";"+data[1]
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
                        # crazy list comprehension in order to identify the name of the atom
                        tmol.create_atom([j[0] for j in pse.items() if j[1][0]==int(line[0])][0],\
                                        float(line[2]),float(line[3]),float(line[4]),'bohr')
                #rest of file has datagrid, x is outer loop, z inner
                tmol.set_vol(nvol,data[6+nat:])
                #finished molecule will be appended to list
                self._mol.append([tmol])

#############################################################################
# WRITE FUNCTIONS
#############################################################################

        def writeFile(self,ftype,mol,filename,param="",coordfmt=""):
                self.outdict[ftype](mol,filename,param,coordfmt)

        def _writeXyz(self,mol,filename,param,coordfmt):
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

        def _writePwi(self,mol,filename,param,coordfmt):
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
                fmt='{0[0][0]:15.10f} {0[0][1]:15.10f} {0[0][2]:15.10f}\n' + \
                    '{0[1][0]:15.10f} {0[1][1]:15.10f} {0[1][2]:15.10f}\n' + \
                    '{0[2][0]:15.10f} {0[2][1]:15.10f} {0[2][2]:15.10f}\n'
                f.write(fmt.format(mol.get_vec()))

                #Close file
                f.close()

class PWParam(dict):

        def __init__(self):
                # make local copy of pse
                self['pse'] = copy.copy(pse)

                # k-point grid
                self['K_POINTS']={'active':'gamma'}
