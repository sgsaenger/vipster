# -*- coding: utf-8 -*-
from . import _indict,_outdict,_paramdict

from copy import deepcopy as _deepcopy
inFormats = _indict.keys()
outFormats = _outdict.keys()

def readFile(filename,fmt):
    """
    Read and parse a given file

    fmt -> file format, needs to be in _indict
    filename -> path to file

    If fmt is in dict, the file will be parsed
    and the resulting trajectory/parameter tuple returned
    """
    with open(filename,'r') as data:
        data=data.readlines()
        return _indict[fmt](filename,data)

def writeFile(mol,fmt,filename,param=None,coordfmt="crystal"):
    """
    Write a given trajectory to disc/stdout

    mol -> molecule to write
    fmt -> file format, needs to be in _outdict
    filename -> path to file
    param -> parameter set to save
    coordfmt -> format in which to save coordinates (bohr/angstrom/crystal/alat)
    """
    if fmt in _paramdict:
        if param is None or param["type"]!=fmt:
            raise Exception(" Writing format "+fmt+" needs accompanying parameter set!")
    with open(filename,'w') as f:
        _outdict[fmt](mol,f,param,coordfmt)

def newParam(prog,var="default"):
    """
    Return a new parameter set for given program

    prog   -> file format, needs to be in _paramdict
    preset -> variant of default set, if present
    """
    return _deepcopy(_paramdict[prog][var])

def availParam(prog=None):
    """
    Return list of available parameter sets

    prog -> if None, return available programs, else return presets
    """
    if prog and prog in _paramdict:
        return _paramdict[prog].keys()
    elif prog:
        return None
    else:
        return _paramdict.keys()
