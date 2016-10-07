# -*- coding: utf-8 -*-
from copy import deepcopy as _deepcopy
from collections import OrderedDict as _ODict
from importlib import import_module as _import
from os import remove
from shutil import move

from vipster.settings import _paramdict

_formats = ["xyz", "pwInput", "pwOutput", "lammpsData", "lammpsCustom", "cube",
            "empire", "aimall", "cpmd", "turbomole", "xsf", "mol2"]
_plugins = []
for i in _formats:
    _plugins.append(_import(".ioplugins." + i, package="vipster"))
_indict = _ODict([(i.argument, i.parser) for i in _plugins])
_outdict = _ODict([(i.argument, i.writer) for i in _plugins if i.writer])
_guiInNames = _ODict([(i.name, i.argument) for i in _plugins])
_guiOutNames = _ODict([(i.name, i.argument) for i in _plugins if i.writer])
_defaultParams = _ODict([(i.argument, i.param) for i in _plugins if i.param])
for i in _defaultParams:
    for j in _defaultParams[i]:
        if j not in _paramdict[i]:
            _paramdict[i][j] = _defaultParams[i][j]
inFormats = _indict.keys()
outFormats = _outdict.keys()


def readFile(filename, fmt):
    """
    Read and parse a given file

    fmt -> file format,  needs to be in _indict
    filename -> path to file

    If fmt is in dict,  the file will be parsed
    and the resulting trajectory/parameter tuple returned
    """
    with open(filename, 'r') as data:
        data = data.readlines()
        return _indict[fmt](filename, data)


def writeFile(mol, fmt, filename, param=None):
    """
    Write a given trajectory to disc/stdout

    mol -> molecule to write
    fmt -> file format,  needs to be in _outdict
    filename -> path to file
    param -> parameter set to save
    """
    if fmt in _paramdict:
        if param is None or param["type"] != fmt:
            raise TypeError(" Writing format " + fmt +
                            " needs accompanying parameter set!")
    # prompt for PermissionError before writing temp file
    t = open(filename, 'w')
    t.close()
    try:
        with open('vipster.tmp', 'w') as f:
            _outdict[fmt](mol, f, param)
    except Exception as e:
        remove('vipster.tmp')
        raise e
    else:
        move('vipster.tmp', filename)


def newParam(prog, var="default"):
    """
    Return a new parameter set for given program

    prog   -> file format,  needs to be in _paramdict
    preset -> variant of default set,  if present
    """
    return _deepcopy(_paramdict[prog][var])


def availParam(prog=None):
    """
    Return list of available parameter sets

    prog -> if None,  return available programs,  else return presets
    """
    if prog and prog in _paramdict:
        return _paramdict[prog].keys()
    elif prog:
        return None
    else:
        return _paramdict.keys()

__all__ = ["inFormats", "outFormats", "readFile", "writeFile",
           "newParam", "availParam"]
