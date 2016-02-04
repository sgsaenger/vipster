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
from . import turbomole

from ..settings import _paramdict

#setup format-lookup-lists
formats=[xyz,pwInput,pwOutput,lammpsData,lammpsCustom,cube,empire,aimall,cpmd,turbomole]
_indict=_ODict([(i.argument,i.parser) for i in formats])
_outdict=_ODict([(i.argument,i.writer) for i in formats if i.writer])
_guiInNames =_ODict([(i.name,i.argument) for i in formats])
_guiOutNames=_ODict([(i.name,i.argument) for i in formats if i.writer])
_defaultParams=_ODict([(i.argument,i.param) for i in formats if i.param])
for i in _defaultParams:
    for j in _defaultParams[i]:
        if j not in _paramdict[i]:
            _paramdict[i][j]=_defaultParams[i][j]

#import members into clean namespace
from . import wrapper
from .wrapper import *
__all__=[i for i in dir(wrapper) if i[0]!='_']
