#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from collections import OrderedDict
from plane_tool import Plane
from vol_tool import Volume
from pick_tool import Picker
from script_tool import Script
from cellmod_tool import CellMod

tools = OrderedDict([
    ('Pick',Picker),
    ('Script',Script),
    ('Cell Mod.',CellMod),
    ('Plane',Plane),
    ('Volume',Volume)])
