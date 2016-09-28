# -*- coding: utf-8 -*-

from collections import OrderedDict

from vipster.gui.tools.plane_tool import Plane
from vipster.gui.tools.vol_tool import Volume
from vipster.gui.tools.pick_tool import Picker
from vipster.gui.tools.script_tool import Script
from vipster.gui.tools.cellmod_tool import CellMod

tools = OrderedDict([
    ('Pick', Picker),
    ('Script', Script),
    ('Cell Mod.', CellMod),
    ('Plane', Plane),
    ('Volume', Volume)])
