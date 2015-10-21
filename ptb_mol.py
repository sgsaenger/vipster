# -*- coding: utf-8 -*-

from math import sqrt
from collections import OrderedDict
from os.path import dirname,expanduser
from json import load, dump

from molecule import Molecule,Trajectory
from ftypeplugins import cli_indict,gui_indict,gui_outdict

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
            self.config = dict()
            self.pse = OrderedDict()
            self.cli_indict = cli_indict
            self.gui_indict = gui_indict
            self.gui_outdict= gui_outdict
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
        with open(dirname(__file__)+'/default.json') as f:
            self.default = load(f,object_pairs_hook=OrderedDict)
        try:
            with open(expanduser('~/.toolbox.json')) as f:
                conf = load(f,object_pairs_hook=OrderedDict)
        except:
            from copy import deepcopy
            conf = deepcopy(self.default)
        self.pse=conf['PSE']
        self.config=conf['General']

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

    def saveConfig(self):
        """Write config and PSE to json-file"""
        with open(expanduser('~/.toolbox.json'),'w') as f:
            dump(OrderedDict([('PSE',self.pse),('General',self.config)]),f,indent=2)
