#include "moltemplate.py.h"
#include <pybind11/eval.h>

py::object Vipster::Plugins::moltemplate(){
py::exec(R"py(
from moltemplate.ltemplify import Ltemplify
import vipster as vp
from io import StringIO

class __vipster_moltemplate:
    name = "moltemplate"
    command = "mol"
    extension = "lt"

    parameters = {
        'name': ("", "Name of the generated template."
                     "If no name is provided, the structure will not be wrapped in a template."
                     "If a forcefield is selected, a name MUST be provided.")
    }

    preset = {
        'atom_style': ((4, ("angle", "atomic", "bond", "charge", "full", "molecular")),
                       "The atom_style to be used by Lammps.\nSee https://lammps.sandia.gov/doc/atom_style.html for details."),
        'forcefield': ((0, ("None", "GAFF", "GAFF2", "OPLSAA", "LOPLSAA", "COMPASS", "TraPPE")),
                       "Prepare a calculation with the given forcefield.")
    }

    @staticmethod
    def writer(mol, p, c, idx):
        # prepare ltemplify with given settings
        args = []
        name = p['name']
        ff = c['forcefield']

        if ff == "None":
            if name:
                args += ['-name', name]
        else:
            if not name:
                raise ValueError("Moltemplate: Use of a forcefield requires providing a template-name.")

            args += ['-ignore-masses', '-ignore-coeffs', '-ignore-bond-types']
            if ff == "GAFF":
                args += ['-name', name+' inherits GAFF',
                         '-preamble', 'import gaff.lt'
                        ]
            elif ff == "GAFF2":
                args += ['-name', name+' inherits GAFF2',
                         '-preamble', 'import gaff2.lt'
                        ]
            elif ff == "OPLSAA":
                args += ['-name', name+' inherits OPLSAA',
                         '-preamble', 'import oplsaa.lt'
                        ]
            elif ff == "LOPLSAA":
                args += ['-name', name+' inherits OPLSAA',
                         '-preamble', 'import loplsaa.lt'
                        ]
            elif ff == "COMPASS":
                args += ['-name', name+' inherits COMPASS',
                         '-preamble', 'import compass_published.lt'
                        ]
            elif ff == "TraPPE":
                args += ['-name', name+' inherits TraPPE',
                         '-preamble', 'import trappe1998.lt'
                        ]

        ltemp = Ltemplify(args)

        # prepare lammps plugin
        lmp = vp.Plugins.lmp
        preset = lmp.makePreset()
        preset['style'] = c['atom_style']
        preset['bonds'] = True

        # create intermediate lammps datafile
        in_file = StringIO(lmp.writer(mol, index=idx, preset=preset))

        # convert via ltemplify
        out_file = StringIO()
        ltemp.Convert(out_file, in_file)

        # return as string
        out_file.seek(0)
        return out_file.read()
)py");
return py::globals()["__vipster_moltemplate"];
}
