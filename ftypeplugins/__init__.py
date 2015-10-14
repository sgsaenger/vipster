from collections import OrderedDict

from . import xyz
from . import pwInput
from . import pwOutput
from . import lammpsData
from . import lammpsCustom
from . import cube
from . import empire
from . import aimall

formats=[xyz,pwInput,pwOutput,lammpsData,lammpsCustom,cube,empire,aimall]
cli_indict=OrderedDict([(i.argument,i.parser) for i in formats])
gui_indict=OrderedDict([(i.name,i.parser) for i in formats])
gui_outdict=OrderedDict([(i.name,i.writer) for i in formats if i.writer])
