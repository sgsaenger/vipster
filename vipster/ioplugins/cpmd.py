# -*- coding: utf-8 -*-
from vipster.molecule import Molecule
from vipster.settings import config

from re import split as regex
from numpy import cos, sqrt
from collections import OrderedDict
from copy import deepcopy

name = 'CPMD Input'
extension = 'cpi'
argument = 'cpmd'

param = {"default": OrderedDict([("type", "cpmd"),
                                 ("name", "New CPMD parameters"),
                                 ("&INFO", ""),
                                 ("&CPMD", ""),
                                 ("&SYSTEM", ""),
                                 ("&PIMD", ""),
                                 ("&PATH", ""),
                                 ("&ATOMS", ""),
                                 ("&DFT", ""),
                                 ("&PROP", ""),
                                 ("&BASIS", ""),
                                 ("&RESP", ""),
                                 ("&PTDDFT", ""),
                                 ("&LINRES", ""),
                                 ("&HARDNESS", ""),
                                 ("&TDDFT", ""),
                                 ("&QMMMM", ""),
                                 ("&CLAS", ""),
                                 ("&EXTE", ""),
                                 ("&VDW", "")])}


def parser(name, data):
    """ Parse CPMD Input file """
    tmol = Molecule(name)
    tparam = deepcopy(param["default"])
    tparam["name"] = name
    i = 0
    ibrav = '14'
    symmetries = {'ISOLATED': '0',
                  'CUBIC': '1',
                  'FACE CENTERED CUBIC': '2',
                  'FCC': '2',
                  'BODY CENTERED CUBIC': '3',
                  'BCC': '3',
                  'HEXAGONAL': '4',
                  'TRIGONAL': '5',
                  'RHOMBOHEDRAL': '5',
                  'TETRAGONAL': '6',
                  'BODY CENTERED TETRAGONAL': '7',
                  'BCT': '7',
                  'ORTHOROMBIC': '8',
                  'MONOCLINIC': '12',
                  'TRICLINIC': '14'}
    fmt = None
    unitvec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    scalevec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    types = []
    while i < len(data):
        if data[i][0] == '!' or data[i].split() == []:
            i += 1
        elif '&' in data[i] and data[i].strip() != "&END":
            # start recording namelist
            nl = data[i].strip()
            i += 1
            while "&END" not in data[i]:
                if nl == "&SYSTEM" and "ANGSTROM" in data[i]:
                    fmt = "angstrom"
                elif nl == "&SYSTEM" and "SCALE" in data[i]:
                    if "CARTESIAN" in data[i]:
                        fmt = "alat"
                    else:
                        fmt = "crystal"
                    for token in data[i].strip().split():
                        if 'S=' in token:
                            scalevec[0][0] = 1 / float(token.strip('S='))
                            scalevec[1][1] = 1 / float(token.strip('S='))
                            scalevec[2][2] = 1 / float(token.strip('S='))
                        elif 'SX=' in token:
                            scalevec[0][0] = 1 / float(token.strip('SX='))
                        elif 'SY=' in token:
                            scalevec[1][1] = 1 / float(token.strip('SY='))
                        elif 'SZ=' in token:
                            scalevec[2][2] = 1 / float(token.strip('SZ='))
                elif nl == "&SYSTEM" and "KPOINTS" in data[i]:
                    if "MONKHORST" in data[i]:
                        kpoints = data[i + 1].split()
                        if "SHIFT=" in data[i]:
                            kpoints += data[i].split("SHIFT=")[1].split()[0:3]
                        else:
                            kpoints += ['0', '0', '0']
                        tmol.setKpoints('active', 'mpg')
                        tmol.setKpoints('mpg', kpoints)
                        i += 1
                    else:
                        kpoints = []
                        crystal = False
                        bands = False
                        if "SCALED" in data[i]:
                            crystal = True
                        if "BANDS" in data[i]:
                            bands = True
                            j = 1
                            line = data[i + j].split()
                            while not all([float(k) == 0 for k in line]):
                                kpoints.append(line[1:4] + [line[0]])
                                kpoints.append(line[4:] + ['0'])
                                j += 1
                                line = data[i + j].split()
                            i += j
                        else:
                            nk = int(data[i + 1])
                            for j in range(nk):
                                kpoints.append(data[i + j + 2].split())
                            i += nk
                        tmol.setKpoints('active', 'discrete')
                        tmol.setKpoints('discrete', kpoints)
                        tmol.setKpoints('options',
                                        {'crystal': crystal, 'bands': bands})
                elif nl == "&SYSTEM" and "SYMMETRY" in data[i]:
                    arg = data[i + 1].strip()
                    if arg.isdigit():
                        ibrav = arg
                    else:
                        ibrav = symmetries[arg]
                    i += 1
                elif nl == "&SYSTEM" and "CELL" in data[i]:
                    cell = list(map(float, data[i + 1].split()))
                    if "VECTORS" in data[i]:
                        while len(cell) < 9:
                            i += 1
                            cell += list(map(float, data[i + 1].split()))
                        ibrav = '-2'
                    else:
                        # cell = [a,b/a,c/a,cos(alpha),cos(beta),cos(gamma)]
                        a, b, c, alpha, beta, gamma = cell
                        tmol.setCellDim(a)
                        if "ABSOLUTE" in data[i]:
                            # cell = [a, b, c, ...]
                            b = b / a
                            c = c / a
                        elif "DEGREE" in data[i]:
                            # cell = [..., alpha, beta, gamma]
                            alpha, beta, gamma =\
                                map(cos, [alpha, beta, gamma])
                    i += 1
                elif nl == "&ATOMS" and data[i][0] == '*':
                    atype = regex('[-_.]', data[i][1:].strip())[0]
                    types.append(atype)
                    tmol.pse[atype]['CPPP'] = data[i][1:].strip()
                    tmol.pse[atype]['CPNL'] = data[i + 1].strip()
                    nat = int(data[i + 2])
                    for j in range(nat):
                        tmol.newAtom(atype, data[i + 3 + j].split()[:3])
                    i += 2 + nat
                elif nl == "&ATOMS" and "CONSTRAINTS" in data[i]:
                    # parse coordinate-fixes
                    j = 1
                    constraints = ""
                    line = data[i + j]
                    while "END" not in line:
                        if "FIX" in line and ("ALL" in line or "QM" in line):
                            for at in range(tmol.nat):
                                tmol.setAtom(at, fix=[True, True, True])
                        elif "FIX" in line and "ELEM" in line:
                            j += 1
                            vals = data[i + j].split()
                            Z = int(vals.pop(0))
                            if "SEQ" in line:
                                beg = int(vals.pop(0)) - 1
                                end = int(vals.pop(0))
                                r = range(beg, end)
                            else:
                                r = range(tmol.nat)
                            for at in r:
                                if tmol.pse[tmol.getAtom(at)[0]]["Z"] == Z:
                                    tmol.setAtom(at, fix=[True, True, True])
                        elif "FIX" in line and "PPTY" in line:
                            j += 1
                            vals = data[i + j].split()
                            PP = int(vals.pop(0)) - 1
                            if "SEQ" in line:
                                beg = int(vals.pop(0)) - 1
                                end = int(vals.pop(0))
                                r = range(beg, end)
                            else:
                                r = range(tmol.nat)
                            for at in r:
                                if tmol.getAtom(at)[0] == types[PP]:
                                    tmol.setAtom(at, fix=[True, True, True])
                        elif "FIX" in line and "SEQ" in line:
                            j += 1
                            vals = data[i + j].split()
                            for at in range(int(vals[0]) - 1, int(vals[1])):
                                tmol.setAtom(at, fix=[True, True, True])
                        elif "FIX" in line and "ATOM" in line:
                            j += 1
                            vals = data[i + j].split()
                            count = int(vals.pop(0))
                            for k in range(count):
                                if not vals:
                                    j += 1
                                    vals = data[i + j].split()
                                tmol.setAtom(int(vals.pop(0)) - 1,
                                             fix=[True, True, True])
                        elif "FIX" in line and "COOR" in line:
                            j += 1
                            vals = data[i + j].split()
                            count = int(vals.pop(0))
                            for k in range(count):
                                j += 1
                                fixcoord = list(map(int, data[i + j].split()))
                                tmol.setAtom(fixcoord[0] - 1, fix=[
                                    bool(not fixcoord[1]),
                                    bool(not fixcoord[2]),
                                    bool(not fixcoord[3])])
                        elif line.strip():
                            constraints += line
                        j += 1
                        line = data[i + j]
                    if constraints:
                        tparam[nl] += "CONSTRAINTS\n" + constraints +\
                                      "END CONSTRAINTS\n"
                    i += j
                elif nl == "&ATOMS" and "ISOTOPE" in data[i]:
                    j = 0
                    for t in types:
                        j += 1
                        mass = data[i + j].strip()
                        tmol.pse[t]['m'] = float(mass)
                    i += j
                else:
                    tparam[nl] += data[i]
                i += 1
        else:
            i += 1
    # apply cell vec (scale SX etc.)
    if scalevec != unitvec:
        tmol.setVec(scalevec, scale=True)
        tmol.setVec(unitvec)
    # parse cell parameters
    if ibrav == '-2':
        tmol.setVec([cell[0:3], cell[3:6], cell[6:9]])
    if ibrav == '1':
        # simple cubic
        pass
    elif ibrav == '2':
        # face centered cubic
        tmol.setVec([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]])
    elif ibrav == '3':
        # body centered cubic
        tmol.setVec([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]])
    elif ibrav == '4':
        # hexagonal
        tmol.setVec([[1, 0, 0], [-0.5, sqrt(3) * 0.5, 0], [0, 0, c]])
    elif ibrav == '5':
        # trigonal
        tx = sqrt((1 - alpha) / 2)
        ty = sqrt((1 - alpha) / 6)
        tz = sqrt((1 + 2 * alpha) / 3)
        tmol.setVec([[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]])
    elif ibrav == '6':
        # simple tetragonal
        tmol.setVec([[1, 0, 0], [0, 1, 0], [0, 0, c]])
    elif ibrav == '7':
        # body centered tetragonal
        tmol.setVec([[0.5, -0.5, c * 0.5],
                     [0.5, 0.5, c * 0.5],
                     [-0.5, -0.5, c * 0.5]])
    elif ibrav == '0' or ibrav == '8':
        # simple orthorhombic
        tmol.setVec([[1, 0, 0], [0, b, 0], [0, 0, c]])
    elif ibrav == '12':
        # simple monoclinic
        tmol.setVec([[1, 0, 0],
                     [b * alpha, b * sqrt(1 - alpha**2), 0],
                     [0, 0, c]])
    elif ibrav == '14':
        # triclinic
        singam = sqrt(1 - gamma**2)
        tmol.setVec([[1, 0, 0], [b * gamma, b * singam, 0],
                     [c * beta, c * (alpha - beta * gamma) / singam,
                      c * sqrt(1 + 2 * alpha * beta * gamma - alpha * alpha -
                               beta * beta - gamma * gamma) / singam]])
    # scale atoms/cell
    if fmt:
        tmol.setFmt(fmt, scale=True)
        tmol.setCellDim(tmol.getCellDim(), fmt=fmt)
    return tmol, tparam


def writer(mol, f, param):
    """
    Save CPMD input file

    Needs both mol and param
    """
    for i in param:
        # always write &CPMD namelist
        if i == "&CPMD":
            f.write("&CPMD\n")
            f.write(param[i])
            f.write("&END\n")
        # expand &SYSTEM with atom and cell information
        elif i == "&SYSTEM":
            f.write("&SYSTEM\n")
            # specify coordfmt
            coordfmt = mol.getFmt()
            if coordfmt == "angstrom":
                f.write("  ANGSTROM\n")
            elif coordfmt == "crystal":
                f.write("  SCALE\n")
            elif coordfmt == "alat":
                f.write("  SCALE CARTESIAN\n")
            # cell vectors are always given explicitely
            f.write("  CELL VECTORS\n")
            vec = mol.getVec() * mol.getCellDim()
            for j in vec:
                f.write("  {:10.5f} {:10.5f} {:10.5f}\n".format(*j))
            # write K-Points if not gamma
            if mol.getKpoints('active') == 'mpg':
                kpoints = mol.getKpoints('mpg')
                f.write("  KPOINTS MONKHORST-PACK")
                if all([float(j) == 0 for j in kpoints[3:]]):
                    f.write("\n")
                else:
                    f.write("  SHIFT={:4s} {:4s} {:4s}\n".format(*kpoints[3:]))
                f.write("  {:4s} {:4s} {:4s}\n".format(*kpoints[:3]))
            elif mol.getKpoints('active') == 'discrete':
                kpoints = mol.getKpoints('discrete')
                opts = mol.getKpoints('options')
                f.write("  KPOINTS")
                if opts["crystal"]:
                    f.write("  SCALED")
                if opts["bands"]:
                    f.write("  BANDS\n")
                    for j in range(0, len(kpoints), 2):
                        line = kpoints[j] + kpoints[j + 1]
                        f.write("  {3:4s} {0:4s} {1:4s}\
                                {2:4s} {4:4s} {5:4s} {6:4s}\n".
                                format(*line))
                    f.write("  0 0. 0. 0. 0. 0. 0.\n")
                else:
                    f.write("\n  " + str(len(kpoints)) + "\n")
                    for j in kpoints:
                        f.write("  {:4s} {:4s} {:4s} {:4s}\n".format(*j))
            # write rest of &SYSTEM namelist
            f.write(param[i])
            f.write("&END\n")
        # put coordinates and PP informations here
        elif i == "&ATOMS":
            f.write("&ATOMS\n")
            atoms = mol.getAtoms(fix=True)
            types = mol.getTypes()
            fixat = []
            fixcoord = []
            count = 0
            for j in types:
                pp = mol.pse[j]['CPPP']
                if not pp:
                    pp = j + config['Default CPMD PP-suffix']
                nl = mol.pse[j]['CPNL']
                if not nl:
                    nl = config['Default CPMD Nonlocality']
                f.write("*" + pp + "\n")
                f.write("  " + nl + "\n")
                f.write("  " + str([a[0] for a in atoms].count(j)) + "\n")
                for k in atoms:
                    if k[0] == j:
                        count += 1
                        f.write("  {:10.5f} {:10.5f} {:10.5f}\n".format(*k[1]))
                        if all(k[2]):
                            fixat.append(count)
                        elif any(k[2]):
                            fixcoord.append((count, k[2]))
            # write rest of &ATOMS namelist
            if not fixat and not fixcoord:
                f.write(param[i])
            else:
                if "CONSTRAINTS" in param[i]:
                    f.write(param[i].partition("CONSTRAINTS")[0])
                f.write("CONSTRAINTS\n")
                if fixat:
                    f.write("FIX ATOMS\n")
                    f.write(str(len(fixat)) + '\n')
                    for j, at in enumerate(fixat):
                        f.write(" {:4d}".format(at))
                        if not (j + 1) % 10:
                            f.write("\n")
                    f.write("\n")
                if fixcoord:
                    f.write("FIX COORDS\n")
                    f.write(str(len(fixcoord)) + '\n')
                    for at in fixcoord:
                        f.write("{:4d}  {:1d} {:1d} {:1d}\n".
                                format(at[0], *at[1]))
                if "CONSTRAINTS" in param[i]:
                    f.write(param[i].partition("CONSTRAINTS")[2])
                else:
                    f.write("END CONSTRAINTS\n")
                    f.write(param[i])
            f.write("ISOTOPE\n")
            for t in types:
                f.write(str(mol.pse[t]['m']) + "\n")
            f.write("&END\n")
        # write other namelists only when they're not empty
        elif i[0] == "&":
            if param[i]:
                f.write(i + "\n")
                f.write(param[i])
                f.write("&END\n")
