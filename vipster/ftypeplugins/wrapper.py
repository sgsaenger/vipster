# -*- coding: utf-8 -*-
from . import _cli_indict,_gui_indict,_cli_outdict,_gui_outdict,_param_dict

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
            return _gui_indict[fmt](filename,data)
        else:
            return _cli_indict[fmt](filename,data)

def writeFile(mol,fmt,filename,param="",coordfmt="",mode="cli"):
    """
    Write a given trajectory to disc/stdout

    mol -> molecule to write
    fmt -> file format, needs to be in outdict
    filename -> path to file
    param -> parameter set to save
    coordfmt -> format in which to save coordinates (bohr/angstrom/crystal/alat)
    """
    with open(filename,'w') as f:
        if mode=="gui":
            _gui_outdict[fmt](mol,f,param,coordfmt)
        else:
            _cli_outdict[fmt](mol,f,param,coordfmt)

def newParam(prog,var="default"):
    """
    Return a new parameter set for given program

    prog   -> file format, needs to be in param_dict
    preset -> variant of default set, if present
    """
    return _param_dict[prog][var]
