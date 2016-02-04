# -*- coding: utf-8 -*-
"""
Parse PSE and config data.
Tries user-config, if not found, falls back to default
"""

from collections import OrderedDict as _ODict
from os.path import dirname as _dirname,expanduser as _expanduser
from json import load as _load,dump as _dump

with open(_dirname(__file__)+'/default.json') as _f:
    default = _load(_f,object_pairs_hook=_ODict)
try:
    with open(_expanduser('~/.vipster.json')) as _f:
        _cfile = _load(_f,object_pairs_hook=_ODict)
except:
    from copy import deepcopy as _deepcopy
    _cfile = _deepcopy(default)
pse = _cfile['PSE']
config = _cfile['General']
_paramdict = _cfile['Parameters']

def saveConfig():
    """Write config and PSE to json-file"""
    with open(_expanduser('~/.vipster.json'),'w') as f:
        _dump(_ODict([('PSE',pse),('General',config),('Parameters',_paramdict)]),f,indent=2)
