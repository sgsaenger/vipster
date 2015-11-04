# -*- coding: utf-8 -*-

from collections import OrderedDict as _ODict
from os.path import dirname as _dirname,expanduser as _expanduser
from json import load as _load,dump as _dump

"""
Parse PSE and config data.
Tries user-config, if not found, falls back to default
"""
with open(_dirname(__file__)+'/default.json') as f:
    default = _load(f,object_pairs_hook=_ODict)
try:
    with open(_expanduser('~/.toolbox.json')) as f:
        cfile = _load(f,object_pairs_hook=_ODict)
except:
    from copy import deepcopy as _deepcopy
    cfile = _deepcopy(default)
pse = cfile['PSE']
config = cfile['General']
del f,cfile

def saveConfig():
    """Write config and PSE to json-file"""
    with open(_expanduser('~/.toolbox.json'),'w') as f:
        _dump(_ODict([('PSE',self.pse),('General',self.config)]),f,indent=2)

"""plugins loaded after config because they might depend on PSE"""
from .ftypeplugins import cli_indict,gui_indict,cli_outdict,gui_outdict,param_dict

def readFile(filename,fmt,mode="cli"):
    """
    Read and parse a given file

    fmt -> file format, needs to be in indict
    filename -> path to file
    mode -> "gui" or "cli", choosing the dictionary

    If fmt is in dict, the file will be parsed
    and the resulting trajectory/parameter tuple returned
    """
    with open(filename,'r') as data:
        data = data.readlines()
        if mode=="gui":
            return gui_indict[fmt](filename,data)
        else:
            return cli_indict[fmt](filename,data)

def writeFile(traj,fmt,filename,param="",coordfmt="",mode="gui"):
    """
    Write a given trajectory to disc/stdout

    traj -> trajectory/molecule to write
    fmt -> file format, needs to be in outdict
    filename -> path to file
    param -> parameter set to save
    coordfmt -> format in which to save coordinates (bohr/angstrom/crystal/alat)
    """
    with open(filename,'w') as f:
        if mode=="gui":
            gui_outdict[fmt](traj,f,param,coordfmt)
        else:
            cli_outdict[fmt](traj,f,param,coordfmt)

def newParam(prog,var="default"):
    """
    Return a new parameter set for given program

    prog   -> file format, needs to be in param_dict
    preset -> variant of default set, if present
    """
    return param_dict[prog][var]
