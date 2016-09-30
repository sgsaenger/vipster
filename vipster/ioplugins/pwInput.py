# -*- coding: utf-8 -*-
from vipster.molecule import Molecule
from vipster.settings import config

from collections import OrderedDict
from numpy import sqrt
from re import split as regex

name = 'PWScf Input'
extension = 'pwi'
argument = 'pwi'

param = {"default": {"type": "pwi",
                     "name": "New PW parameters",
                     "&control": OrderedDict([('calculation', 'scf')]),
                     "&system": OrderedDict(),
                     "&electrons": OrderedDict()}}


def parser(name, data):
    """ Parse PWScf input files

    Namelists will be parsed and saved in PW parameter set
    Supported Cards:
    - ATOMIC_SPECIES
    - ATOMIC_POSITIONS
    - K_POINTS
    - CELL_PARAMETERS
    Not supported:
    - CONSTRAINTS
    - OCCUPATIONS
    - ATOMIC_FORCES (PWSCFv5)
    """
    tmol = Molecule(name)
    tparam = {"type": "pwi", "name": name}
    tvec = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    # parse data and create tparam
    while data:
        header = data.pop(0).strip().split()
        #  ignore empty lines
        if not header:
            pass
        #  parse namelists
        elif header[0][0] == '&':
            tnl = OrderedDict()
            #  parse entries
            line = regex(', \s*(?![^()]*\))', data.pop(0).strip())
            while line[0] != '/':
                for j in range(len(line)):
                    if line[j]:
                        tnl[line[j].split('=')[0].strip()] =\
                            line[j].split('=')[1].strip()
                line = regex(', \s*(?![^()]*\))', data.pop(0).strip())
            tparam[header[0].lower()] = tnl
        #  parse card
        elif header[0][0].isupper():
            # ATOMIC_SPECIES:
            # Name   Weight  PP-file
            if header[0] == 'ATOMIC_SPECIES':
                for i in range(int(tparam['&system']['ntyp'])):
                    line = data.pop(0).strip().split()
                    tmol.pse[line[0]]['m'] = float(line[1])
                    tmol.pse[line[0]]['PWPP'] = line[2]
            # ATOMIC_POSITIONS fmt
            # Name   x   y   z
            elif header[0] == 'ATOMIC_POSITIONS':
                fmt = header[1].strip('{()}')
                tmol.newAtoms(int(tparam['&system']['nat']))
                for i in range(int(tparam['&system']['nat'])):
                    # support empty lines
                    temp = data.pop(0).split()
                    while not temp:
                            temp = data.pop(0).split()
                    tmol.setAtom(i, temp[0], temp[1:4],
                                 fix=[not int(x) for x in temp[4:]])
            # K_POINTS fmt
            elif header[0] == 'K_POINTS':
                active = header[1].strip('{()}')
                # Gamma point only
                if active == 'gamma':
                    tmol.setKpoints('active', 'gamma')
                # Monkhorst Pack Grid:
                # x y z offset
                elif active == 'automatic':
                    line = data.pop(0).split()
                    tmol.setKpoints('mpg', line)
                    tmol.setKpoints('active', 'mpg')
                # else:
                # number of kpoints
                # x y z weight
                else:
                    nk = int(data.pop(0).split()[0])
                    kpoints = []
                    for i in range(nk):
                        kpoints.append(data.pop(0).split())
                    tmol.setKpoints('discrete', kpoints)
                    tmol.setKpoints('active', 'discrete')
                    if active == 'tpiba':
                        pass
                    elif active == 'tpiba_b':
                        tmol.setKpoints('options', {'bands': True,
                                                    'crystal': False})
                    elif active == 'crystal':
                        tmol.setKpoints('options', {'bands': False,
                                                    'crystal': True})
                    elif active == 'crystal_b':
                        tmol.setKpoints('options', {'bands': True,
                                                    'crystal': True})
            # CELL_PARAMETERS
            # only needed if ibrav=0
            # tbd changed between pw4 and pw5, ignored for now
            elif header[0] == 'CELL_PARAMETERS':
                for i in [0, 1, 2]:
                    line = data.pop(0).strip().split()
                    tvec[i] = [float(x) for x in line]
            else:
                pass
    # Identify cell parameter representation, parse it
    tmol.setCellDim(tparam['&system']['celldm(1)'])
    ibrav = tparam['&system']['ibrav']
    if ibrav == '0':
        # check if CELL_PARAMETERS card has been read
        # if not present, throw error
        if tvec == [[0, 0, 0], [0, 0, 0], [0, 0, 0]]:
                raise ValueError('ibrav=0, but CELL_PARAMETERS missing')
        else:
                tmol.setVec(tvec)
    elif ibrav == '1':
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
        c = float(tparam['&system']['celldm(3)'])
        tmol.setVec([[1, 0, 0], [-0.5, sqrt(3) * 0.5, 0], [0, 0, c]])
    elif ibrav == '5':
        # trigonal
        alpha = float(tparam['&system']['celldm(4)'])
        tx = sqrt((1 - alpha) / 2)
        ty = sqrt((1 - alpha) / 6)
        tz = sqrt((1 + 2 * alpha) / 3)
        tmol.setVec([[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]])
    elif ibrav == '-5':
        # trigonal, alternative
        alpha = float(tparam['&system']['celldm(4)'])
        tx = sqrt((1 - alpha) / 2)
        ty = sqrt((1 - alpha) / 6)
        tz = sqrt((1 + 2 * alpha) / 3)
        u = (tz - 2 * sqrt(2) * ty) / sqrt(3)
        v = (tz + sqrt(2) * ty) / sqrt(3)
        tmol.setVec([[u, v, v], [v, u, v], [v, v, u]])
    elif ibrav == '6':
        # simple tetragonal
        c = float(tparam['&system']['celldm(3)'])
        tmol.setVec([[1, 0, 0], [0, 1, 0], [0, 0, c]])
    elif ibrav == '7':
        # body centered tetragonal
        c = float(tparam['&system']['celldm(3)'])
        tmol.setVec([[0.5, -0.5, c * 0.5],
                     [0.5, 0.5, c * 0.5],
                     [-0.5, -0.5, c * 0.5]])
    elif ibrav == '8':
        # simple orthorhombic
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        tmol.setVec([[1, 0, 0], [0, b, 0], [0, 0, c]])
    elif ibrav == '9':
        # basis centered orthorhombic
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        tmol.setVec([[0.5, b * 0.5, 0], [-0.5, b * 0.5, 0], [0, 0, c]])
    elif ibrav == '10':
        # face centered orthorhombic
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        tmol.setVec([[0.5, 0, c * 0.5],
                     [0.5, b * 0.5, 0],
                     [0, b * 0.5, c * 0.5]])
    elif ibrav == '11':
        # body centered orthorhombic
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        tmol.setVec([[0.5, b * 0.5, c * 0.5],
                     [-0.5, b * 0.5, c * 0.5],
                     [-0.5, -b * 0.5, c * 0.5]])
    elif ibrav == '12':
        # simple monoclinic
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        gamma = float(tparam['&system']['celldm(4)'])
        tmol.setVec([[1, 0, 0],
                     [b * gamma, b * sqrt(1 - gamma**2), 0],
                     [0, 0, c]])
    elif ibrav == '-12':
        # simple monoclinic, alternate definition
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        beta = float(tparam['&system']['celldm(5)'])
        tmol.setVec([[1, 0, 0], [0, b, 0],
                     [c * beta, 0, c * sqrt(1 - beta**2)]])
    elif ibrav == '13':
        # base centered monoclinic
        c = float(tparam['&system']['celldm(3)'])
        b = float(tparam['&system']['celldm(2)'])
        gamma = float(tparam['&system']['celldm(4)'])
        tmol.setVec([[0.5, 0, -c * 0.5],
                     [b * gamma, b * sqrt(1 - gamma**2), 0],
                     [0.5, 0, c * 0.5]])
    elif ibrav == '14':
        # triclinic
        b = float(tparam['&system']['celldm(2)'])
        c = float(tparam['&system']['celldm(3)'])
        gamma = float(tparam['&system']['celldm(6)'])
        beta = float(tparam['&system']['celldm(5)'])
        alpha = float(tparam['&system']['celldm(4)'])
        singam = sqrt(1 - gamma**2)
        tmol.setVec([[1, 0, 0], [b * gamma, b * singam, 0],
                    [c * beta, c * (alpha - beta * gamma) / singam,
                     c * sqrt(1 + 2 * alpha * beta * gamma - alpha * alpha -
                              beta * beta - gamma * gamma) / singam]])
    # create atoms after creating cell:
    tmol.setFmt(fmt, scale=True)
    # delete nat, ntype and celldm before returning
    del tparam['&system']['nat']
    del tparam['&system']['ntyp']
    del tparam['&system']['celldm(1)']
    return tmol, tparam


def writer(mol, f, param):
    """
    Save PWScf input file

    Needs both mol and param
    """
    # &control, &system and &electron namelists are mandatory
    for i in ['&control', '&system', '&electrons']:
        f.write(i + '\n')
        # write all parameters
        if i == '&system':
            f.write(' nat=' + str(mol.nat) + '\n')
            f.write(' ntyp=' + str(mol.ntyp) + '\n')
            f.write(' celldm(1)=' + str(mol.getCellDim(fmt='bohr')) + '\n')
        for j in range(len(param[i])):
            f.write(' ' + list(param[i].keys())[j] + '=' +
                    list(param[i].values())[j] + '\n')
        f.write('/\n\n')
    # &ions only when needed
    if '&ions' in param:
        f.write('&ions' + '\n')
        for j in range(len(param['&ions'])):
            f.write(' ' + list(param['&ions'].keys())[j] + '=' +
                    list(param['&ions'].values())[j] + '\n')
        f.write('/\n\n')
    elif param['&control']['calculation'] in\
            ["'relax'", "'vc-relax'", "'md'", "'vc-md'"]:
        raise KeyError('&ions namelist required, but not present')
    # &cell only when needed
    if '&cell' in param:
        f.write('&cell' + '\n')
        for j in range(len(param['&cell'])):
            f.write(' ' + list(param['&cell'].keys())[j] + '=' +
                    list(param['&cell'].values())[j] + '\n')
        f.write('/\n\n')
    elif param['&control']['calculation'] in ["'vc-relax'", "'vc-md'"]:
        raise KeyError('&cell namelist required, but not present')
    # ATOMIC_SPECIES card:
    f.write('ATOMIC_SPECIES' + '\n')
    types = list(mol.getTypes())
    for i in range(len(mol.getTypes())):
        atom = types[i]
        pp = str(mol.pse[atom]['PWPP'])
        if not pp:
            pp = atom + config['Default PWScf PP-suffix']
        f.write(atom + '    ' + str(mol.pse[atom]['m']) + '   ' + pp + '\n')
    f.write('\n')
    # ATOMIC_POSITIONS
    f.write('ATOMIC_POSITIONS' + ' ' + mol.getFmt() + '\n')
    for i in range(mol.nat):
        atom = mol.getAtom(i, fix=True)
        if any(atom[-1]):
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f} {:1d} {:1d} {:1d}'.
                    format(atom[0], atom[1][0], atom[1][1], atom[1][2],
                           not atom[-1][0], not atom[-1][1],
                           not atom[-1][2]) + '\n')
        else:
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                atom[0], atom[1][0], atom[1][1], atom[1][2]) + '\n')
    f.write('\n')
    # K_POINTS
    active = mol.getKpoints('active')
    # Gamma point only
    if active == 'gamma':
        f.write('K_POINTS gamma\n')
    # MPGrid:
    # x y z offset
    elif active == 'mpg':
        f.write('K_POINTS automatic\n')
        auto = mol.getKpoints('mpg')
        f.write('{:4s} {:4s} {:4s} {:4s} {:4s} {:4s}'.format(
                auto[0], auto[1], auto[2], auto[3], auto[4], auto[5]) + '\n')
    # number of kpoints
    # x y z weight
    else:
        opts = mol.getKpoints('options')
        if opts['bands']:
            if opts['crystal']:
                active = 'crystal_b'
            else:
                active = 'tpiba_b'
        elif opts['crystal']:
            active = 'crystal'
        else:
            active = 'tpiba'
        f.write('K_POINTS ' + active + '\n')
        disc = mol.getKpoints('discrete')
        f.write(str(len(disc)) + '\n')
        for i in range(len(disc)):
            f.write('{:4s} {:4s} {:4s} {:4s}'.format(
                    disc[i][0], disc[i][1], disc[i][2], disc[i][3]) + '\n')
    f.write('\n')
    # Cell parameters
    f.write('CELL_PARAMETERS' + '\n')
    fmt = '{0[0][0]:15.10f} {0[0][1]:15.10f} {0[0][2]:15.10f}\n' + \
          '{0[1][0]:15.10f} {0[1][1]:15.10f} {0[1][2]:15.10f}\n' + \
          '{0[2][0]:15.10f} {0[2][1]:15.10f} {0[2][2]:15.10f}\n'
    f.write(fmt.format(mol.getVec()))