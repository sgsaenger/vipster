# -*- coding: utf-8 -*-
from vipster.molecule import Molecule
from collections import OrderedDict

name = "Lammps Data"
extension = 'lmp'
argument = 'lmp'
param = {"default": {"type": "lmp",
                     "name": "New LAMMPS",
                     "atom_style": "atomic",
                     "bonds": False, "angles": False,
                     "dihedrals": False, "impropers": False}}
# Mapping:
# 0 = atom-ID (int)
# 1 = molecule-ID (int)
# 2 = atom-type (int)
# 3,4,5 = x,y,z (int)
# 6 = charge (float)
lammps_atom_style = OrderedDict([
    ("angle", "{0:d} {1:d} {2:d} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("atomic", "{0:d} {2:d} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("bond", "{0:d} {1:d} {2:d} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("charge", "{0:d} {2:d} {6:s} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("full", "{0:d} {1:d} {2:d} {6:s} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("molecular", "{0:d} {1:d} {2:d} {3: .4f} {4: .4f} {5: .4f}\n")])
"""unsupported due to unsupported properties:
    ("body", "{0:d} {2:d} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("dipole", "{0:d} {2:d} {6:s} {3: .4f} {4: .4f} {5: .4f} 0 0 0\n"),
    ("dpd", "{0:d} {2:d} {} {3: .4f} {4: .4f} {5: .4f}\n")
    ("electron", "{0:d} {2:d} {6:s} 0 0 {3: .4f} {4: .4f} {5: .4f}\n"),
    ("ellipsoid", "{0:d} {2:d} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("line", "{0:d} {1:d} {2:d} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("meso", "{0:d} {2:d} {} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("peri", "{0:d} {2:d} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("smd", "{0:d} {2:d} {1:d} {} {} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("sphere", "{0:d} {2:d} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("template", "{0:d} {1:d} {} {} {2:d} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("tri", "{0:d} {1:d} {2:d} {} {} {3: .4f} {4: .4f} {5: .4f}\n"),
    ("wavepacket", "{0:d} {2:d} {} {} {} {} {} {} {3: .4f} {4: .4f} {5: .4f}\n") # noqa: E501
"""


def parser(name, data):
    """ Parse Lammps data file

    Element parsed via comment in 'Masses' part
    Assumes angstrom, atom_style in [angle,atomic,bond,charge,full,molecular]
    """
    tmol = Molecule(name)
    i = 0
    tvec = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    while i < len(data):
        line = data[i].strip()
        if 'atoms' in line:
            nat = int(line.split()[0])
        elif 'atom types' in line:
            types = [0] * int(line.split()[0])
        elif 'xlo xhi' in line:
            tvec[0][0] = float(line.split()[1]) - float(line.split()[0])
        elif 'ylo yhi' in line:
            tvec[1][1] = float(line.split()[1]) - float(line.split()[0])
        elif 'zlo zhi' in line:
            tvec[2][2] = float(line.split()[1]) - float(line.split()[0])
        elif 'xy xz yz' in line:
            tvec[1][0], tvec[2][0], tvec[2][1] = line.split()[:3]
        elif 'Masses' in line:
            for j in range(i + 2, i + 2 + len(types)):
                if '#' in data[j]:
                    el = data[j].split('#')
                    t = el[1].strip()
                    types[int(el[0].split()[0]) - 1] = t
                    tmol.pse[t]['m'] = float(el[0].split()[1])
                else:
                    raise NotImplementedError("Can't deduce elements")
            i += len(types) + 1
        elif 'Atoms' in line:
            atomlist = []
            for j in range(i + 2, i + 2 + nat):
                atomlist.append(data[j].strip().split())
            if '#' in line and line.split('#')[1].strip() in lammps_atom_style:
                style = lammps_atom_style[line.split('#')[1].
                                          strip()].split('} ')
                atypepos = style.index('{2:d')
                cpos = style.index('{3: .4f')
                try:
                    qpos = style.index('{6:s')
                except:
                    qpos = None
            else:
                # make a guess for column-content:
                nargs = len(atomlist[0])
                atypepos = 1
                cpos = -3
                testimageflags = 0
                for j in range(nargs, 1, -1):
                    try:
                        test = [int(k[j]) for k in atomlist]
                    except:
                        pass
                    else:
                        # check for imageflags (3 integer columns)
                        if nargs > 7 and j > 3:
                            testimageflags += 1
                            if testimageflags == 3:
                                cpos = -6
                        # check if one of the first two entries
                        # and not an empty charge column
                        elif j == 2 and set(test) != set([0]):
                            atypepos = j
                            break
                # guess where charge is found:
                if nargs + cpos - 1 > atypepos:
                    qpos = atypepos + 1
                else:
                    qpos = False
            # gotta load 'em all!
            tmol.newAtoms(nat)
            for j, at in enumerate(atomlist):
                tmol.setAtom(j, types[int(at[atypepos]) - 1],
                             at[cpos:cpos + 3 if cpos + 3 else None],
                             at[qpos] if qpos else None)
            tmol.setFmt('angstrom', scale=True)
            tmol.setVec(tvec)
            tmol.setCellDim(1, fmt='angstrom')
        i += 1
    return tmol, None


def writer(mol, f, param):
    """
    Save Lammps input data

    Cell parameters have to be in lower triangular form
    Assumes angstrom
    """
    # calculate necessary values:
    if param["bonds"] or param["angles"]\
            or param["dihedrals"] or param["impropers"]\
            or param["atom_style"] in ["angle", "bond", "full", "line",
                                       "molecular", "smd", "template", "tri"]:
        bondlist = []
        for i in mol.getBonds():
            for j in i:
                bondlist.append(tuple(sorted(j[0:2])))
        bondlist = sorted(set(bondlist))
        bondtypes = []
        for i in bondlist:
            bondtypes.append(tuple(sorted([mol.getAtom(j)[0] for j in i])))
        bondtypelist = sorted(list(set(bondtypes)))
    if param["angles"] or param["dihedrals"] or param["impropers"]:
        anglelist = []
        for i in range(len(bondlist)):
            for j in range(i + 1, len(bondlist)):
                if bondlist[i][0] == bondlist[j][0]:
                    anglelist.append(tuple([bondlist[i][1],
                                            bondlist[i][0],
                                            bondlist[j][1]]))
                elif bondlist[i][1] == bondlist[j][0]:
                    anglelist.append(tuple([bondlist[i][0],
                                            bondlist[i][1],
                                            bondlist[j][1]]))
                elif bondlist[i][1] == bondlist[j][1]:
                    anglelist.append(tuple([bondlist[i][0],
                                            bondlist[i][1],
                                            bondlist[j][0]]))
        angletypes = []
        for i in anglelist:
            at = [mol.getAtom(j)[0] for j in i]
            if at[0] > at[2]:
                at.reverse()
            angletypes.append(tuple(at))
        angletypelist = sorted(list(set(angletypes)))
    if param["dihedrals"]:
        dihedrallist = []
        for i in range(len(anglelist)):
            for j in range(i + 1, len(anglelist)):
                a1 = anglelist[i]
                a2 = anglelist[j]
                if a1[1] == a2[0] and a2[1] == a1[2]:
                    dihedrallist.append(tuple([a1[0], a1[1], a1[2], a2[2]]))
                elif a1[1] == a2[0] and a2[1] == a1[0]:
                    dihedrallist.append(tuple([a1[2], a1[1], a1[0], a2[2]]))
                elif a1[1] == a2[2] and a2[1] == a1[2]:
                    dihedrallist.append(tuple([a1[0], a1[1], a1[2], a2[0]]))
                elif a1[1] == a2[2] and a2[1] == a1[0]:
                    dihedrallist.append(tuple([a1[2], a1[1], a1[0], a2[0]]))
        dihedraltypes = []
        for i in dihedrallist:
            dt = [mol.getAtom(j)[0] for j in i]
            if dt[0] > dt[3]:
                dt.reverse()
            elif dt[0] == dt[3] and dt[1] > dt[2]:
                dt.reverse()
            dihedraltypes.append(tuple(dt))
        dihedtypelist = sorted(list(set(dihedraltypes)))
    if param["impropers"]:
        from itertools import combinations
        improperlist = []
        for i in range(mol.nat):
            neigh = []
            for j in bondlist:
                if i in j:
                    neigh.append(j[not j.index(i)])
            for j in combinations(neigh, 3):
                improperlist.append((i,) + j)
        impropertypes = []
        for i in improperlist:
            it = [mol.getAtom(j)[0] for j in i]
            impropertypes.append(tuple(it))
        imptypelist = sorted(list(set(impropertypes)))
    moleculeid = [0] * mol.nat
    if param["atom_style"] in ["angle", "bond", "full", "line",
                               "molecular", "smd", "template", "tri"]:
        molbondlist = [set(i) for i in bondlist]

        def groupsets(setlist):
            tlist = [setlist[0]]
            for i in setlist[1:]:
                matched = False
                for j in tlist:
                    if i & j:
                        j.update(i)
                        matched = True
                if not matched:
                    tlist.append(i)
            if tlist == setlist:
                return tlist
            else:
                return groupsets(tlist)
        moleculelist = groupsets(molbondlist)
        for i, m in enumerate(moleculelist):
            for at in m:
                moleculeid[at] = i
    # write header:
    f.write('\n' + str(mol.nat) + ' atoms\n')
    f.write(str(mol.ntyp) + ' atom types\n')
    if param["bonds"] and bondlist:
        f.write(str(len(bondlist)) + ' bonds\n')
        f.write(str(len(bondtypelist)) + ' bond types\n')
        for i, j in enumerate(bondtypelist):
            f.write('#{:d} {:s} {:s}\n'.format(i + 1, *j))
    if param["angles"] and anglelist:
        f.write(str(len(anglelist)) + ' angles\n')
        f.write(str(len(angletypelist)) + ' angle types\n')
        for i, j in enumerate(angletypelist):
            f.write('#{:d} {:s} {:s} {:s}\n'.format(i + 1, *j))
    if param["dihedrals"] and dihedrallist:
        f.write(str(len(dihedrallist)) + ' dihedrals\n')
        f.write(str(len(dihedtypelist)) + ' dihedral types\n')
        for i, j in enumerate(dihedtypelist):
            f.write('#{:d} {:s} {:s} {:s} {:s}\n'.format(i + 1, *j))
    if param["impropers"] and improperlist:
        f.write(str(len(improperlist)) + ' impropers\n')
        f.write(str(len(imptypelist)) + ' improper types\n')
        for i, j in enumerate(imptypelist):
            f.write('#{:d} {:s} {:s} {:s} {:s}\n'.format(i + 1, *j))
    f.write('\n')
    # if cell is orthogonal,  write vectors:
    v = mol.getVec() * mol.getCellDim(fmt='angstrom')
    if not v.diagonal(1).any() and not v.diagonal(2).any():
        f.write('{: .4f} {: .4f} xlo xhi\n'.format(0.0, v[0][0]))
        f.write('{: .4f} {: .4f} ylo yhi\n'.format(0.0, v[1][1]))
        f.write('{: .4f} {: .4f} zlo zhi\n'.format(0.0, v[2][2]))
        if v.diagonal(-1).any() or v.diagonal(-2).any():
            f.write('{: .4f} {: .4f} {: .4f} xy xz yz\n'.format(v[1][0],
                                                                v[2][0],
                                                                v[2][1]))
        f.write('\n')
    else:
        raise ValueError('Cell not in suitable Lammps format')
    # Masses section: (always)
    f.write('Masses\n\n')
    atomtypes = list(mol.getTypes())
    for i, j in enumerate(atomtypes):
        f.write('{:d} {:2.4f} #{:s}\n'.format(i + 1, mol.pse[j]['m'], j))
    # Atoms section: (always)
    f.write('\nAtoms # {:}\n\n'.format(param["atom_style"]))
    string = lammps_atom_style[param["atom_style"]]
    for i, at in enumerate(mol.getAtoms(charge=True, fmt='angstrom')):
        f.write(string.format(i + 1, moleculeid[i] + 1,
                              atomtypes.index(at[0]) + 1,
                              at[1][0], at[1][1], at[1][2], at[-1]))

    def convertIndices(x):
        return list(map(lambda y: tuple(map(lambda z: z + 1, y)), x))
    # Bonds section:
    if param["bonds"] and bondlist:
        f.write('\nBonds\n\n')
        bondlist = convertIndices(bondlist)
        for i, j in enumerate(bondlist):
            f.write('{:d} {:d} {:d} {:d}\n'.
                    format(i + 1, bondtypelist.index(bondtypes[i]) + 1, *j))
    # Angles section:
    if param["angles"] and anglelist:
        f.write('\nAngles\n\n')
        anglelist = convertIndices(anglelist)
        for i, j in enumerate(anglelist):
            f.write('{:d} {:d} {:d} {:d} {:d}\n'.
                    format(i + 1, angletypelist.index(angletypes[i]) + 1, *j))
    # Dihedrals section:
    if param["dihedrals"] and dihedrallist:
        f.write('\nDihedrals\n\n')
        dihedrallist = convertIndices(dihedrallist)
        for i, j in enumerate(dihedrallist):
            f.write('{:d} {:d} {:d} {:d} {:d} {:d}\n'.
                    format(i + 1,
                           dihedtypelist.index(dihedraltypes[i]) + 1, *j))
    # Impropers section:
    if param["impropers"] and improperlist:
        f.write('\nImpropers\n\n')
        improperlist = convertIndices(improperlist)
        for i, j in enumerate(improperlist):
            f.write('{:d} {:d} {:d} {:d} {:d} {:d}\n'.
                    format(i + 1, imptypelist.index(impropertypes[i]) + 1, *j))
    f.write('\n')
