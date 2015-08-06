from collections import OrderedDict

import xyz
import pwInput
import pwOutput
import lammpsData
import lammpsCustom
import cube

formats=[xyz,pwInput,pwOutput,lammpsData,lammpsCustom,cube]
cli_indict=OrderedDict([(i.argument,i.parser) for i in formats])
gui_indict=OrderedDict([(i.name,i.parser) for i in formats])
gui_outdict=OrderedDict([(i.name,i.writer) for i in formats if i.writer])
