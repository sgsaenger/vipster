# -*- coding: utf-8 -*-

#read local config and setup global variables/dictionaries
from vipster.settings import *
from vipster.settings import __all__

#main data-class
from vipster.molecule import Molecule
__all__ += ["Molecule"]

#i/o-routines
from vipster.iowrapper import *
from vipster.iowrapper import __all__ as _aio_
__all__ += _aio_

#gui-launcher
from vipster.gui.main import launchVipster
__all__ += ["launchVipster"]
