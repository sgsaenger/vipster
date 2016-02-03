# -*- coding: utf-8 -*-
from collections import OrderedDict as _ODict
from . import xyz
from . import pwInput
from . import pwOutput
from . import lammpsData
from . import lammpsCustom
from . import cube
from . import empire
from . import aimall
from . import cpmd

from ..settings import _userParams

#setup format-lookup-lists
formats=[xyz,pwInput,pwOutput,lammpsData,lammpsCustom,cube,empire,aimall,cpmd]
_indict=_ODict([(i.argument,i.parser) for i in formats])
_outdict=_ODict([(i.argument,i.writer) for i in formats if i.writer])
_guiInNames =_ODict([(i.name,i.argument) for i in formats])
_guiOutNames=_ODict([(i.name,i.argument) for i in formats if i.writer])
_paramdict=_ODict([(i.argument,i.param) for i in formats if i.param])
for i in _userParams:
    for j in _userParams[i]:
        _paramdict[i][j]=_userParams[i][j]

#import members into clean namespace
from . import wrapper
from .wrapper import *
__all__=[i for i in dir(wrapper) if i[0]!='_']
