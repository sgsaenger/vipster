# -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy

from vipster.mol_c import setBondsC, makeVolGradient
from vipster.settings import pse, config

# constants
bohrrad = 0.52917721092
invbohr = 1.88972612565
fmts = ['angstrom', 'bohr', 'crystal', 'alat']
kpfmts = ['gamma', 'mpg', 'discrete']

# per-step properties
_properties = {'_atom_name': [],
               '_atom_coord': np.array([], 'f'),
               '_atom_fix': [],
               '_atom_hidden': [],
               '_atom_charge': [],
               '_bonds': [[], [], [], [], [], [], [], []],
               '_bond_cutfac': None,
               '_bonds_outdated': True,
               '_default_fmt': 'bohr',
               '_selection': [],
               '_script_definitions': dict(),
               '_celldm': 1.0,
               '_vec': np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'f'),
               '_vecinv': np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'f'),
               '_vol': None,
               '_vol_grad': None,
               '_vol_off': None,
               '_iso_val': (False, 0, False),
               '_vol_plane': (False, 0, 0),
               '_mil_plane': (False, (0, 0, 0), (0, 0, 0)),
               'comment': '',
               '_undo_depth': 0,
               '_undo_stack': [],
               '_undo_temp': []}


class _step(object):
    """
    Actual step data

    Includes:
    Atom symbols/coordinates (_atom_name/_atom_coord)
    Bonds between atoms (generated when requested) (_bonds)
    List of fixed atoms in PWScf calculation (_atom_fix)
    XYZ comment line (comment)
    Cube-style volume data (_vol)
    Cell geometry (_celldm/_vec)
    """

    def __init__(self, pse=None):
        self.pse = pse
        for i in _properties:
            setattr(self, i, deepcopy(_properties[i]))
        self._bond_cutfac = float(config['Bond cutoff factor'])

    ######################################################
    # ATOM FUNCTIONS
    ######################################################

    def newAtom(self, name='C', coord=[0., 0., 0.], charge=0.,
                fix=[False, False, False], hidden=False, fmt=None):
        """Make new atom
        name -> element symbol
        coord -> coordinates
        fix -> 0 for no relaxation in PW
        fmt -> str in ['bohr','angstrom','crystal','alat']
        """
        self._atom_name.append(name)
        self._atom_coord.resize([self.nat + 1, 3])
        self._atom_coord[-1] = self._coordToBohr(coord, fmt)
        while len(fix) < 3:
            fix += [False]
        self._atom_fix.append(fix)
        self._atom_hidden.append(hidden)
        self._atom_charge.append(charge)
        self.delSelection()
        self._bonds_outdated = True

    def newAtoms(self, count):
        """Make 'count' new empty atoms"""
        self._atom_name.extend(['X'] * count)
        self._atom_coord.resize([self.nat + count, 3])
        self._atom_fix.extend([[False, False, False]] * count)
        self._atom_hidden.extend([False] * count)
        self._atom_charge.extend([0.] * count)
        self.delSelection()
        self._bonds_outdated = True

    def delAtom(self, index):
        """Remove atom"""
        del self._atom_name[index]
        self._atom_coord = np.delete(self._atom_coord, index, 0)
        del self._atom_fix[index]
        del self._atom_hidden[index]
        del self._atom_charge[index]
        self.delSelection()
        self._bonds_outdated = True

    def setAtom(self, index, name=None, coord=None, charge=None, fix=None,
                hidden=None, fmt=None):
        """Modify a given atom

        index -> atom id
        rest -> new properties
        """
        if name:
            self._atom_name[index] = name
        if coord is not None:
            self._atom_coord[index] = self._coordToBohr(coord, fmt)
        if fix:
            while len(fix) < 3:
                fix += [False]
            self._atom_fix[index] = fix
        if hidden is not None:
            self._atom_hidden[index] = hidden
        if charge is not None:
            self._atom_charge[index] = charge
        self._bonds_outdated = True

    nat = property(lambda self: len(self._atom_coord))

    def getAtoms(self, charge=False, fix=False, hidden=False, fmt=None):
        """Return names and coordinates for all atoms

        fmt -> str in ['angstrom','bohr','crystal','alat']
        """
        l = [self._atom_name, map(lambda f: self._coordFromBohr(f, fmt),
             self._atom_coord)]
        if charge:
            l.append(self._atom_charge)
        if fix:
            l.append(self._atom_fix)
        if hidden:
            l.append(self._atom_hidden)
        return list(zip(*l))

    def getAtom(self, index, charge=False, fix=False, hidden=False, fmt=None):
        """Return exact copy of one atom

        index -> index of atom
        fmt -> str in ['angstrom','bohr','crystal','alat']
        """
        l = [self._atom_name[index],
             self._coordFromBohr(self._atom_coord[index],
             fmt)]
        if charge:
            l.append(self._atom_charge[index])
        if fix:
            l.append(self._atom_fix[index])
        if hidden:
            l.append(self._atom_hidden[index])
        return l

    def getTypes(self):
        """Return types of atoms

        sorted by atomic number
        """
        return sorted(sorted(set(self._atom_name)),
                      key=lambda n: self.pse[n]["Z"])

    ntyp = property(lambda self: len(set(self._atom_name)))

    #######################################################
    # COORD FMT FUNCTIONS
    ######################################################

    def setFmt(self, fmt, scale=False):
        """Set default format for future atom-operations

        fmt -> str in ['bohr','angstrom','crystal','alat']
        scale -> if true, existing atoms are rescaled from fmt to bohr"""
        self._default_fmt = fmt
        if scale:
            if fmt == 'angstrom':
                self._atom_coord *= invbohr
            elif fmt == 'crystal':
                self._atom_coord = np.dot(self._atom_coord, self._vec) * \
                    self._celldm
            elif fmt == 'alat':
                self._atom_coord = self._atom_coord * self._celldm

    def getFmt(self):
        """Return default format"""
        return self._default_fmt

    def _coordToBohr(self, coord, fmt=None):
        """Transform given coordinates to bohr

        coord -> new coordinates in given format
        fmt -> str in ['angstrom','bohr','crystal','alat']
        retval -> new coordinates in bohr
        """
        coord = np.array(coord, 'f')
        if not fmt:
            fmt = self._default_fmt
        if fmt == 'angstrom':
            return coord * invbohr
        elif fmt == 'bohr':
            return coord
        elif fmt == 'crystal':
            return np.dot(coord, self._vec) * self._celldm
        elif fmt == 'alat':
            return coord * self._celldm

    def _coordFromBohr(self, coord, fmt):
        """Transform given coordinates from bohr

        coord -> old coordinates in bohr
        fmt -> str in ['angstrom','bohr','crystal','alat']
        retval -> new coordinates in given format
        """
        if not fmt:
            fmt = self._default_fmt
        if fmt == 'angstrom':
            return coord * bohrrad
        elif fmt == 'bohr':
            return coord
        elif fmt == 'crystal':
            return np.dot(coord, self._vecinv) / self._celldm
        elif fmt == 'alat':
            return coord / self._celldm

    ######################################################
    # CELL FUNCTIONS
    ######################################################

    def setCellDim(self, cdm, scale=False, fmt=None):
        """Set new cell-dimension

        cdm -> new cell-dimension
        scale -> if True: scale coordinates of atoms (fixed crystal-coord.)
        fmt -> str in ['angstrom','bohr']
        """
        if not fmt:
            fmt = self._default_fmt
        if fmt == 'angstrom':
            cdm *= invbohr
        if scale:
            ratio = cdm / self._celldm
            self._atom_coord *= ratio
        self._celldm = float(cdm)
        self._bonds_outdated = True

    def setVec(self, vec, scale=False):
        """Set new cell-vectors

        vec -> 3x3 list of new vectors
        scale -> if True: scale coordinates of atoms (fixed crystal-coord.)
        """
        vec = np.array(vec, 'f')
        inv = self._vecinv
        if scale:
            self._atom_coord = np.dot(np.dot(self._atom_coord, inv), vec)
        self._vec = vec
        self._vecinv = np.linalg.inv(self._vec)
        self._bonds_outdated = True

    def getCellDim(self, fmt=None):
        """Return cell-dimension

        fmt -> str in ['angstrom','bohr','crystal','alat']
        """
        if not fmt:
            fmt = self._default_fmt
        if fmt == 'angstrom':
            return self._celldm * bohrrad
        else:
            return self._celldm

    def getVec(self):
        """Return cell-vectors"""
        return self._vec

    def getCenter(self, com=False):
        """Return center-coordinates of molecule

        if com is True: return center of mass
        else returns center of cell
        """
        if com and self.nat:
            return (np.max(self._atom_coord, axis=0) +
                    np.min(self._atom_coord, axis=0)) / 2
        else:
            return (self._vec[0] + self._vec[1] + self._vec[2]) *\
                self._celldm / 2

    ######################################################
    # UNDO FUNCTIONS
    ######################################################

    def initUndo(self):
        """Save atoms temporarily

        Will be passed to undo-stack if needed
        """
        if not self._undo_depth:
            self._undo_temp = [
                deepcopy(self._atom_name),
                deepcopy(self._atom_coord),
                deepcopy(self._atom_charge),
                deepcopy(self._atom_fix),
                deepcopy(self._atom_hidden),
                deepcopy(self._celldm),
                deepcopy(self._vec),
                deepcopy(self._vecinv)]
        self._undo_depth += 1

    def saveUndo(self, name=''):
        """Put temp-save on stack

        name -> name of action to reverse
        """
        self._undo_depth -= 1
        if not self._undo_depth:
            self._undo_stack.append([name] + self._undo_temp)

    def getUndo(self):
        """Return name of last undo, None otherwise"""
        if self._undo_stack:
            return self._undo_stack[-1][0]
        else:
            return None

    def undo(self):
        """Undo last action

        Pop atom infos from stack if present
        """
        if not self._undo_stack:
            return
        _, self._atom_name, self._atom_coord, self._atom_charge,\
            self._atom_fix, self._atom_hidden, self._celldm,\
            self._vec, self._vecinv = self._undo_stack.pop()
        self.delSelection()
        self._bonds_outdated = True

    ######################################################
    # SELECTION FUNCTIONS
    ######################################################

    def addSelection(self, item):
        """Add an item to selection

        item -> [index,(offset in multiples of vec)]
        """
        if item in self._selection:
            self._selection.remove(item)
        else:
            self._selection.append(item)

    def delSelection(self):
        """Clear selection"""
        self._selection = []

    def getSelection(self):
        """Return the selection"""
        return self._selection

    ######################################################
    # BOND FUNCTIONS
    ######################################################

    def getBonds(self, cutfac=None):
        """Return bonds

        Sets bonds if not present.
        Returns list of lists.
        First entry: Bonds inside of cell
        Other entries: Periodic bonds (see setBonds)
        """
        if not cutfac:
            cutfac = self._bond_cutfac
        if self._bonds_outdated:
            self.setBonds(cutfac)
        elif self._bond_cutfac != cutfac:
            self.setBonds(cutfac)
        return self._bonds

    def setBonds(self, cutfac):
        """Set bonds

        Cutoff criterium:
        Sum of covalent radii times cutfac
        Generates list of lists, each lists contains:
        [idx_at1,idx_at2,offset_at1,offset_at2,distance]
        Actual setting performed by C-routine setBondsC
        """
        self._bonds = [[], [], [], [], [], [], [], []]
        if self.nat < 2:
            return
        at_c = self._atom_coord
        cutoff = np.array([self.pse[i]['bondcut'] for i in self._atom_name],
                          'f')
        n = np.zeros(3, 'f')
        v = self.getVec() * self.getCellDim('bohr')
        off = [[(n, n)],                           # orig
               [(v[0], n)],                        # x
               [(v[1], n)],                        # y
               [(v[0] + v[1], n), (v[0], v[1])],     # xy, x-y
               [(v[2], n)],                        # z
               [(v[0] + v[2], n), (v[0], v[2])],     # xz, x-z
               [(v[1] + v[2], n), (v[1], v[2])],     # yz, y-z
               [(v[0] + v[1] + v[2], n), (v[0] + v[1], v[2]),
               (v[0] + v[2], v[1]), (v[1] + v[2], v[0])]]  # xyz,xy-z,x-yz,-xyz
        self._bonds = setBondsC(at_c, cutoff, cutfac, off)
        self._bond_cutfac = cutfac
        self._bonds_outdated = False

    #####################################################
    # EDIT FUNCTIONS
    #####################################################

    def mult(self, x, y, z):
        """Multiply cell

        x,y,z -> Integer multipliers along corresponding vector
        """
        self.initUndo()
        vec = self._vec * self._celldm
        mult = [x, y, z]
        nat = self.nat
        self._atom_coord.resize([self.nat * np.prod(mult), 3])
        for k in range(3):
            self._atom_name.extend((mult[k] - 1) * self._atom_name)
            self._atom_fix.extend((mult[k] - 1) * self._atom_fix)
            self._atom_hidden.extend((mult[k] - 1) * self._atom_hidden)
            self._atom_charge.extend((mult[k] - 1) * self._atom_charge)
            for i in range(1, mult[k]):
                for j in range(nat):
                    self._atom_coord[j + i * nat] =\
                        self._atom_coord[j] + i * vec[k]
            nat = nat * mult[k]
        self.setVec(self._vec * [[x], [y], [z]])
        self._bonds_outdated = True
        self.saveUndo('multiply cell')

    def crop(self):
        """Crop all atoms outside of cell"""
        self.initUndo()
        nat = self.nat
        dellist = []
        for i in range(nat):
            at = self.getAtom(i, fmt='crystal')
            if np.any(at[1] > 1.000001) or np.any(at[1] < -0.000001) or\
                    np.any(np.isclose(at[1] - 1, 0, atol=1.e-6)):
                dellist.append(i)
        dellist.reverse()
        for i in dellist:
            self.delAtom(i)
        self._bonds_outdated = True
        self.saveUndo('crop atoms')

    def wrap(self):
        """Wrap all atoms outside of cell to the inside"""
        self.initUndo()
        nat = self.nat
        for i in range(nat):
            at = self.getAtom(i, fmt='crystal')
            self.setAtom(i, at[0], at[1] % 1, fmt='crystal')
        self._bonds_outdated = True
        self.saveUndo('wrap atoms')

    def reshape(self, newvec):
        """Reshape cell

        newvec -> new cell vectors

        Transforms cell so new cell is a cut of
        bulk-material with given geometry.
        PBC must be manually conserved!
        """
        self.initUndo()
        self.wrap()
        newdim = abs(np.array(newvec)).sum(0)
        olddim = abs(self._vec).sum(0)
        m = 1
        while np.any(olddim * m < newdim):
            m += 1
        self.mult(m, m, m)
        oldcenter = self.getCenter()
        self.setVec(newvec)
        newcenter = self.getCenter()
        for i in range(self.nat):
            self._atom_coord[i] += newcenter - oldcenter
        self.crop()
        self.saveUndo('reshape cell')

    def align(self, vec, direc):
        """Align cell vectors

        vec -> int specifying the vector to align
        direc -> str specifying the direction to align to
        """

        if direc == 'x':
            d = np.array([1, 0, 0], 'f')
        elif direc == 'y':
            d = np.array([0, 1, 0], 'f')
        elif direc == 'z':
            d = np.array([0, 0, 1], 'f')
        else:
            raise ValueError('Align vectors: invalid direction!')

        v = self.getVec()[int(vec)]
        v = v / np.linalg.norm(v)

        if np.all(np.equal(abs(v), d)):
            return

        self.initUndo()
        ax = np.cross(v, d)
        ax = ax / np.linalg.norm(ax)
        c = np.float(np.dot(v, d))
        theta = np.arccos(np.dot(v, d))
        s = np.float(-np.sin(theta))
        ic = np.float(1. - c)
        mat = np.array(
            [[ic * ax[0] * ax[0] + c, ic * ax[0] * ax[1] - s * ax[2],
              ic * ax[0] * ax[2] + s * ax[1]],
             [ic * ax[0] * ax[1] + s * ax[2], ic * ax[1] * ax[1] + c,
              ic * ax[1] * ax[2] - s * ax[0]],
             [ic * ax[0] * ax[2] - s * ax[1], ic * ax[1] * ax[2] + s * ax[0],
              ic * ax[2] * ax[2] + c]], 'f')
        mat = np.array([np.dot(i, mat) for i in self.getVec()])
        self.setVec(mat, scale=True)
        self._bonds_outdated = True
        self.saveUndo('align cell')

    #####################################################
    # VOLUME DATA FUNCTIONS
    #####################################################

    def setVol(self, vol, off):
        """Set volume data

        vol -> data-grid
        off -> offset of data
        """
        self._vol = np.array(vol, 'f')
        self._vol_off = np.array(off, 'f')
        self._vol_grad = makeVolGradient(self._vol)

    def getVol(self):
        """Return volume data if present"""
        return self._vol

    def getVolOffset(self):
        """Return offset of volume data"""
        return self._vol_off

    def getVolGradient(self):
        """Return volume gradient"""
        return self._vol_grad

    def setIsoVal(self, enabled, val, pm):
        """Set current isovalue"""
        self._iso_val = (enabled, val, pm)

    def getIsoVal(self):
        """Return current isovalue"""
        return self._iso_val

    def setVolPlane(self, enabled, direction, position):
        """Set position/direction of volume-heatmap"""
        self._vol_plane = (enabled, direction, position)

    def getVolPlane(self):
        """Return current volume-heatmap settings"""
        return self._vol_plane

    def setMilPlane(self, enabled, indices):
        """Set miller-indices"""
        self._mil_plane = (enabled, indices)

    def getMilPlane(self):
        """Return miller-indices"""
        return self._mil_plane

    #####################################################
    # SCRIPTING FUNCTIONS
    #####################################################

    def evalScript(self, script):
        """Evaluate script

        script -> list of str containing script

        Operators are defined in ops-dictionary.
        Defined groups of atoms are saved in instanced list '_script_group'.
        Other evaluated operations are placed on stack and
        executed after successful parsing.
        """
        def evalOperator(op, args):
            """
            op -> Given operator
            args -> required arguments for operator
            retval -> closure for evaluation of given arguments
            """
            def getToken(arglist, sep=' '):
                arglist = arglist.split(sep, 1)
                arg = arglist[0]
                if sep != ' ':
                    arg += sep
                if len(arglist) > 1:
                    rest = arglist[1].strip()
                else:
                    rest = ''
                return arg, rest

            def evalOpTokens(arglist):
                res = [op]
                for t in args:
                    # list of atoms
                    if t == 'l':
                        # list of indices
                        if arglist[0] == '[':
                            arg, arglist = getToken(arglist, ']')
                            arg = arg.strip('[]').split(', ')
                            l = []
                            for j in arg:
                                # interpret range
                                if '-' in j:
                                    low, high = j.split('-')
                                    l.extend(range(int(low), int(high) + 1))
                                # atom index
                                else:
                                    l.append(int(j))
                            res.append(l)
                        else:
                            arg, arglist = getToken(arglist)
                            # reference to definition
                            if arg.isalpha():
                                if arg == 'all':
                                    res.append(range(self.nat))
                                elif arg == 'sel':
                                    res.append(set([a[0] for a in
                                               self._selection]))
                                else:
                                    res.append(self._script_definitions[arg])
                            # single index
                            elif arg.isdigit():
                                res.append([int(arg)])
                            elif '-' in arg:
                                low, high = arg.split('-')
                                res.append(range(int(low), int(high) + 1))
                    # float (factor or angle)
                    elif t in 'fF':
                        # nothing left and optional:
                        if not arglist and t == 'F':
                            continue
                        arg, arglist = getToken(arglist)
                        try:
                            arg = float(arg)
                            res.append(float(arg))
                        except Exception as e:
                            if t == 'f':
                                raise ValueError('Invalid float: ' + str(arg))
                    # index for loops
                    elif t == 'i':
                        arg, arglist = getToken(arglist)
                        res.append(int(arg))
                    # arbitrary names for defined groups
                    elif t in 'sS':
                        arg, arglist = getToken(arglist)
                        if t.isupper() and not arg:
                            continue
                        res.append(arg)
                    # valid vector for operations
                    elif t in 'vV':
                        # nothing left and optional:
                        if not arglist and t == 'V':
                            continue
                        # explicit vector
                        if arglist[0] == '(':
                            arg, arglist = getToken(arglist, ')')
                            arg = eval(arg)
                            if len(arg) == 4 and type(arg[3]) is str:
                                arg = self._coordToBohr(arg[0:3], arg[3])
                            elif len(arg) == 3:
                                arg = self._coordToBohr(arg)
                            else:
                                raise ValueError('Invalid vector: ' + str(arg))
                        # implicit vector
                        elif arglist[0].isdigit():
                            arg, arglist = getToken(arglist)
                            # difference between atoms
                            if '-' in arg:
                                arg = arg.split('-')
                                arg = self._atom_coord[int(arg[0])] -\
                                    self._atom_coord[int(arg[1])]
                            # position of atom
                            else:
                                arg = np.array(self._atom_coord[int(arg)])
                        # implicit negative vector
                        elif arglist[0] == '-':
                            arg, arglist = getToken(arglist)
                            arg = np.array(-self._atom_coord[int(arg[1:])])
                        # fail when not a vector and vector needed
                        elif t == 'v':
                            raise ValueError('Not a valid vector: ' + str(arg))
                        # continue when not a vector and not needed
                        else:
                            continue
                        res.append(arg)
                return res
            return evalOpTokens
        # dictionary returns closure containing target operation
        # and argument-check
        # upper-case: optional
        ops = {'rot': evalOperator(self.rotate, 'lfvV'),
               'shi': evalOperator(self.shift, 'lvF'),
               'def': evalOperator('define', 'slSSS'),
               'sel': evalOperator('select', 'lSSS'),
               'mir': evalOperator(self.mirror, 'lvvV'),
               'psh': evalOperator(self.pshift, 'lvl'),
               'par': evalOperator(self.parallelize, 'lll')}
        stack = []
        script = script.strip().split('\n')
        # check for errors,  parse and prepare
        for line in script:
            try:
                operator, args = line.split(' ', 1)
                action = ops[operator[0:3].lower()](args)
            except KeyError as e:
                return 'Wrong Operator: ' + operator
            except Exception as e:
                return e.message
            else:
                if action[0] == 'define':
                    self._script_definitions[action[1]] = \
                        self.filterAtoms(action[2], *action[3:])
                elif action[0] == 'select':
                    self.delSelection()
                    sel = self.filterAtoms(action[1], *action[2:])
                    sel = zip(sel, [(0, 0, 0)] * len(sel))
                    for i in sel:
                        self.addSelection(i)
                else:
                    stack.append(action)
        # delete previous definitions
        self._script_definitions = {}
        # if everything went well,  execute operations
        try:
            self.initUndo()
            for action in stack:
                action[0](*action[1:])
        except Exception as e:
            self._undo_depth = 0
            return e.message
        else:
            self.saveUndo('script')
            return 'Success!'

    def filterAtoms(self, atoms, *filter):
        atoms = set(atoms)
        for f in filter:
            if f[0] in ['x', 'y', 'z']:
                axdir = {'x': 0, 'y': 1, 'z': 2, '<': -1, '>': 1}
                axis = axdir[f[0]]
                direction = axdir[f[1]]
                target = float(f[2:])
                remove = []
                for at in atoms:
                    if direction * (self.getAtom(at)[1][axis] - target) <= 0:
                        remove.append(at)
                for at in remove:
                    atoms.remove(at)
            elif f[0] == 'c':
                target = int(f[1:])
                bonds = self.getBonds()
                remove = []
                for at in atoms:
                    count = 0
                    for b in bonds:
                        for p in b:
                            count += p.count(at)
                    if count != target:
                        remove.append(at)
                for at in remove:
                    atoms.remove(at)
            elif f[0] == 't':
                remove = []
                for at in atoms:
                    if self.getAtom(at)[0] != f[1:]:
                        remove.append(at)
                for at in remove:
                    atoms.remove(at)
        return atoms

    def rotate(self, atoms, angle, ax, shift=np.zeros(3, 'f')):
        """Rotate group of atoms

        atoms -> list of atoms
        angle -> angle in degree
        ax -> rotation-axis
        shift -> optional vector shifting the rotation-axis
        """
        self.initUndo()
        angle = np.radians(angle)
        c = np.float(np.cos(angle))
        s = np.float(-np.sin(angle))
        ic = np.float(1. - c)
        if np.isclose(np.linalg.norm(ax), 0, atol=1.e-6):
            raise ValueError("0 Vector can't be rotation axis")
        ax = ax / np.linalg.norm(ax)
        mat = np.array(
            [[ic * ax[0] * ax[0] + c, ic * ax[0] * ax[1] - s * ax[2],
              ic * ax[0] * ax[2] + s * ax[1]],
             [ic * ax[0] * ax[1] + s * ax[2], ic * ax[1] * ax[1] + c,
              ic * ax[1] * ax[2] - s * ax[0]],
             [ic * ax[0] * ax[2] - s * ax[1], ic * ax[1] * ax[2] + s * ax[0],
              ic * ax[2] * ax[2] + c]], 'f')
        for i in atoms:
            self._atom_coord[i] =\
                np.dot(self._atom_coord[i] - shift, mat) + shift
        self._bonds_outdated = True
        self.saveUndo('rotate atoms')

    def shift(self, atoms, vector, factor=1.0):
        """Shift group of atoms

        atoms -> list of atoms
        vector -> shift-vector
        """
        self.initUndo()
        for i in atoms:
            self._atom_coord[i] += np.array(vector, 'f') * factor
        self._bonds_outdated = True
        self.saveUndo('shift atoms')

    def mirror(self, atoms, v1, v2, shift=np.zeros(3, 'f')):
        """Mirror group of atoms

        atoms -> list of atoms
        v1,v2 -> vectors defining the plane
        shift -> optional vector shifting the mirror-plane
        """
        self.initUndo()
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)
        for i in atoms:
            pos = self._atom_coord[i]
            proj = np.dot(pos - shift, normal) * normal
            self._atom_coord[i] = pos - 2 * proj
        self._bonds_outdated = True
        self.saveUndo('mirror atoms')

    def parallelize(self, atoms, p1, p2):
        """Parallelize planes

        atoms -> list of atoms
        p1 -> plane in molecule/mobile target
              second element is static point of reference
        p2 -> plane on surface/static target
        """
        self.initUndo()
        if len(p1) != 3:
            raise ValueError('Mobile plane not valid')
        if len(p2) != 3:
            raise ValueError('Static plane not valid')
        p1 = [self._atom_coord[i] for i in p1]
        p2 = [self._atom_coord[i] for i in p2]
        d11 = p1[0] - p1[1]
        d12 = p1[2] - p1[1]
        n1 = np.cross(d11, d12)
        n1 = n1 / np.linalg.norm(n1)
        d21 = p2[0] - p2[1]
        d22 = p2[2] - p2[1]
        n2 = np.cross(d21, d22)
        n2 = n2 / np.linalg.norm(n2)
        angle = np.degrees(np.arccos(np.dot(n2, n1)))
        if np.isclose(angle, 0, atol=1.e-6) \
                or np.isclose(angle, 0, atol=1.e-6):
            return
        axis = np.cross(n1, n2)
        self.rotate(atoms, angle, axis, p1[1])
        self.saveUndo('parallelize planes')

    def pshift(self, atoms, vec, plane):
        """Shift over plane

        atoms -> list of atoms
        vec -> shift vector
              x/y parallel to static plane
              z perpendicular
        pl -> plane on surface/static target
              pl[0]-pl[1] give x-coordinate
        """
        self.initUndo()
        if len(plane) != 3:
            raise ValueError('Static plane not valid')
        plane = [self._atom_coord[i] for i in plane]
        x = plane[0] - plane[1]
        t = plane[2] - plane[1]
        z = np.cross(x, t)
        z = z / np.linalg.norm(z)
        shift = np.dot(vec, [x, np.cross(x, z), z])
        self.shift(atoms, shift)
        self.saveUndo('shift over plane')


class Molecule(_step):
    """
    Main data container

    consistent interface for single molecules
    and trajectories

    Includes:
    Interface to modify all cells at once
    K-Point settings (_kpoints)
    PSE-Overlay over central settings (pse)
    """
    def __init__(self, name="New Mol", steps=1):
        """
        Create container with `steps` empty molecules
        """
        self.name = name
        self.pse = self._pse(pse)
        self.steps = [_step(self.pse) for i in range(steps)]
        if steps:
            self.curStep = 0
        else:
            self.curStep = None
        self._kpoints = {'active': 'gamma',
                         'mpg': ['1', '1', '1', '0', '0', '0'],
                         'discrete': [],
                         'options': {'crystal': False, 'bands': False}}

    def copy(self):
        """
        Return copy of current step
        """
        m = Molecule(name="Copy of " + self.name, steps=0)
        m.steps.append(deepcopy(self.steps[self.curStep]))
        m.curStep = 0
        return m

    def copyAll(self):
        """
        Return copy of trajectory
        """
        m = deepcopy(self)
        m.name = "Copy of " + self.name
        return m

    def __len__(self):
        return len(self.steps)

    def newStep(self):
        """
        Create new step
        """
        self.steps.append(_step(self.pse))
        self.changeStep(len(self.steps) - 1)

    def copyStep(self, num=None):
        """
        Copy given or current step
        """
        if num is None:
            num = self.curStep
        self.steps.append(deepcopy(self.steps[num]))
        self.changeStep(len(self.steps) - 1)

    def changeStep(self, num):
        """
        Select step 'num'
        """
        self.curStep = num

    def setVecAll(self, vec, scale=False):
        """
        Set new cell-vectors for all included steps
        """
        if scale:
            for step in self.steps:
                step.setVec(vec, scale)
        else:
            vec = np.array(vec, 'f')
            vecinv = np.linalg.inv(self._vec)
            for step in self.steps:
                step._vec = vec
                step._vecinv = vecinv
                step._bonds_outdated = True

    def setCellDimAll(self, cdm, scale=False, fmt=None):
        """
        Set new cell-dimension for all included steps
        """
        for step in self.steps:
            step.setCellDim(cdm, scale, fmt)

    def setKpoints(self, mode, kpoints):
        """Set and modify active k-points

        mode -> str in [active,automatic,disc]
        kpoints -> corresponding argument
        self._kpoints['active'] in ['gamma','mpg','discrete']
        self._kpoints['mpg'] = [x,y,z,xoff,yoff,zoff]
        self._kpoints['discrete'] = [[x,y,z,w],...]
        self._kpoints['options'] = {'crystal':False,'bands':False}
        """
        self._kpoints[mode] = kpoints

    def getKpoints(self, mode):
        """Return active k-points or k-point-settings"""
        return self._kpoints[mode]

    ######################################################
    # LOCAL PSE-DICT
    # overlay for global dict
    ######################################################

    class _pse(dict):
        """Special dictionary

        Upon requesting information for an atom-type,
        determines likely referenced type and makes
        local copy
        """
        def __init__(self, cpse):
            super(Molecule._pse, self).__init__()
            self.cpse = cpse

        def __getitem__(self, key):
            if key not in self:
                if key in self.cpse:
                    self[key] = deepcopy(self.cpse[key])
                    source = key
                else:
                    for i in range(len(key), 0, -1):
                        source = None
                        if key[:i] in self.cpse:
                            source = key[:i]
                        elif key[:i].upper() in self.cpse:
                            source = key[:i].upper()
                        if source:
                            self[key] = deepcopy(self.cpse[source])
                            break
                    if key not in self:
                        self[key] = deepcopy(self.cpse['X'])
                        source = 'X'
                if not self[key]["PWPP"]:
                    self[key]["PWPP"] = source +\
                        config["Default PWScf PP-suffix"]
                if not self[key]["CPPP"]:
                    self[key]["CPPP"] = source +\
                        config["Default CPMD PP-suffix"]
            return super(Molecule._pse, self).__getitem__(key)

# add step-properties to container:
for i in _properties:
    def getter(self, i=i):
        return self.steps[self.curStep].__getattribute__(i)

    def setter(self, val, i=i):
        self.steps[self.curStep].__setattr__(i, val)

    setattr(Molecule, i, property(getter, setter))
