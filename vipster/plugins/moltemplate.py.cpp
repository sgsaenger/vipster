#include "moltemplate.py.h"
#include <pybind11/eval.h>

py::object Vipster::Plugins::moltemplate(){
py::exec(R"py(
from moltemplate.ltemplify import Ltemplify
from moltemplate.lttree import LttreeSettings, LttreeParseArgs,\
    StaticObj, InstanceObj, BasicUI, ExecCommands, VarRef
import vipster as vp
from io import StringIO
from collections import defaultdict

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
    def parser(name, file):
        # prepare settings and initial (re-)parsing of file
        settings = LttreeSettings()
        LttreeParseArgs(['', name], settings, main=True)

        # prepare tree for templates and their instances
        g_objectdefs = StaticObj('', None)
        g_objects = InstanceObj('', None)

        # the command buffers for moltemplate
        g_static_commands = []
        g_instance_commands = []

        # parse file(s) into command buffers and object trees
        BasicUI(settings,
                g_objectdefs, g_objects,
                g_static_commands, g_instance_commands)

        def mkLmpMol(name, atoms, types, cell=None, bonds=None, btypes=''):
            lmp = vp.Plugins.lmp
            natom = len([at for at in atoms.split('\n') if at.split('#')[0].strip()])
            ntype = len([at for at in types.split('\n') if at.split('#')[0].strip()])
            if not cell:
                disableCell = True
                cell = '0 500000 xlo xhi\n'\
                       '0 500000 ylo yhi\n'\
                       '0 500000 zlo zhi\n'
            else:
                disableCell = False
            if not bonds:
                tmpfile = f"""
{natom} atoms
{ntype} atom types

{cell}

Masses

{types}

Atoms

{atoms}
"""
            else:
                nbond = len([b for b in bonds.split('\n') if b.split('#')[0].strip()])
                nbtyp = len([b for b in btypes.split('\n')])
                tmpfile = f"""
{natom} atoms
{ntype} atom types
{nbond} bonds
{nbtyp} bond types
{btypes}

{cell}

Masses

{types}

Atoms

{atoms}

Bonds

{bonds}
"""
            mol = lmp.parser(name, tmpfile)
            mol[0].comment = name
            if disableCell:
                mol[0].enableCell(False)
            return mol

        if 'atom' not in g_objects.categories:
            # template-only file
            if not g_objectdefs.children:
                raise ValueError("Moltemplate: neither templates nor actual atoms present in file, aborting")
            mol = vp.Molecule(name, 0)
            # recursively convert templates into steps
            def parseTemplate(mol, t_name, t_impl):
                for t in t_impl.children.items():
                    parseTemplate(mol, *t)
                if t_impl.commands == [[]]:
                    return
                cell = None
                types = ''
                atoms = ''
                natom = 0
                for cmd in t_impl.commands[:-1]:
                    if cmd.filename == 'Data Boundary':
                        cell = cmd.tmpl_list[0].text
                    elif cmd.filename == 'Data Masses':
                        refs = []
                        def cmdToString(c):
                            nonlocal refs
                            if type(c) is VarRef:
                                refs.append(c)
                                return c.binding.value
                            else:
                                return c.text
                        types = ''.join(map(cmdToString, cmd.tmpl_list))
                        # if types are well-behaved, create our own type-annotations
                        tmp = [t for t in types.split('\n') if t.split('#')[0].strip()]
                        if len(tmp) == len(refs):
                            types = ''.join(['{} # {}\n'.format(t, n.nptr.leaf_node.name) for t,n in zip(tmp, refs)])
                for cmd in t_impl.commands[-1]:
                    if cmd.filename == 'Data Atoms':
                        def cmdToString(c):
                            nonlocal natom
                            if type(c) is VarRef:
                                if type(c) is VarRef:
                                    if c.descr_str[:4] == "atom":
                                        if c.prefix == "@":
                                            return c.binding.value
                                        else:
                                            natom += 1
                                            return str(natom)
                                    else:
                                        return '0'
                                else:
                                    return ''
                            else:
                                return c.text
                        atoms = ''.join(map(cmdToString, cmd.tmpl_list))
                    elif cmd.filename == 'Data Bonds':
                        # TODO: not really possible atm?
                        # can only identify atoms via string representation
                        # which may break if equivalent but mismatching names are used
                        pass

                # parse template via intermediate lammps file
                mol.newStep(mkLmpMol(t_name, atoms, types, cell)[0])
            for t in g_objectdefs.children.items():
                parseTemplate(mol, *t)
            return mol
        else:
            # file instantiates templates, simulate full lttree-execution
            tmp = defaultdict(list)
            ExecCommands(g_static_commands, tmp, settings, True)
            ExecCommands(g_instance_commands, tmp, settings, True)

            # convert sections to single strings
            types = ''.join(tmp['Data Masses'])
            atoms = ''.join(tmp['Data Atoms'])
            bonds = ''.join(tmp['Data Bonds'])
            cell = ''.join(tmp['Data Boundary'])
            btypes = ''.join(['# {} {}\n'.format(i, b.nptr.leaf_node.name) for i,b in enumerate(g_objectdefs.categories['bond'].bindings.values(), 1)])

            return mkLmpMol(name, atoms, types, cell, bonds, btypes)

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
