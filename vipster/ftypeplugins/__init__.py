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
from . import turbomole

#setup format-lookup-lists
formats=[xyz,pwInput,pwOutput,lammpsData,lammpsCustom,cube,empire,aimall,turbomole]
_cli_indict=_ODict([(i.argument,i.parser) for i in formats])
_cli_outdict=_ODict([(i.argument,i.writer) for i in formats if i.writer])
_gui_indict=_ODict([(i.name,i.parser) for i in formats])
_gui_outdict=_ODict([(i.name,i.writer) for i in formats if i.writer])
_param_dict=_ODict([(i.argument,i.param) for i in formats if i.param])

#import members into clean namespace
from . import wrapper
from .wrapper import *
__all__=[i for i in dir(wrapper) if i[0]!='_']
