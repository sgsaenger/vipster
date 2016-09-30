# -*- coding: utf-8 -*-
from vipster.molecule import Molecule

name = 'xsf'
extension = 'xsf'
argument = 'xsf'
param = None
writer = None


def parser(name, data):
    """ Parse (animated) XCrysDen structure files

    Optionally contains datagrids/animations/forces
    Parses multiple datagrids into multiple steps
    """
    animated = False
    periodic = False
    pbcformats = ["MOLECULE", "POLYMER", "SLAB", "CRYSTAL"]
    i = 0
    while i < len(data):
        if data[i][0] == '#':
            pass
        elif "ANIMSTEPS" in data[i]:
            tmol = Molecule(name, steps=int(data[i].strip("ANIMSTEPS \r\n")))
            animated = True
        elif data[i].strip() in pbcformats:
            periodic = True
        # Parse coordinates of isolated molecule
        elif not periodic and "ATOMS" in data[i]:
            if animated:
                tmol.changeStep(int(data[i].strip("ATOMS \r\n")) - 1)
            else:
                tmol = Molecule(name, steps=1)
            j = 1
            while i + j < len(data) and len(data[i + j].split()) in [4, 7]:
                line = data[i + j].split()
                tmol.newAtom(line[0], line[1:4])
                j += 1
            tmol.setFmt("angstrom", scale=True)
            i += j - 1
        # Parse periodic structure containing cell and coordinates
        elif periodic and "PRIMVEC" in data[i]:
            vec = [data[i + j].split() for j in range(1, 4)]
            if animated and data[i].strip("PRIMVEC \r\n"):
                tmol.changeStep(int(data[i].strip("PRIMVEC ")) - 1)
                tmol.setVec(vec)
                tmol.setCellDim(1, fmt='angstrom')
            elif animated:
                tmol.setVecAll(vec)
                tmol.setCellDimAll(1, fmt='angstrom')
            else:
                tmol = Molecule(name, steps=1)
                tmol.setVec(vec)
                tmol.setCellDim(1, fmt='angstrom')
            i += 3
        elif periodic and "PRIMCOORD" in data[i]:
            if animated:
                tmol.changeStep(int(data[i].strip("PRIMCOORD \r\n")) - 1)
            nat = int(data[i + 1].split()[0])
            tmol.newAtoms(nat)
            for j in range(nat):
                line = data[j + 2 + i].split()
                tmol.setAtom(j, line[0], line[1:4])
            tmol.setFmt("angstrom", scale=True)
            i += nat + 1
        elif periodic and "CONVVEC" in data[i]:
            # could be used, but not atm
            pass
        elif periodic and "CONVCOORD" in data[i]:  # pragma: no cover
            # should not even be there
            pass
        elif "BEGIN_" in data[i]:
            datatype = data[i].strip("BEGIN_BLOCK_ \r\n")
            endblock = "END_BLOCK_" + datatype
            endgrid = "END_" + datatype
            j = 1
            grids = []
            # only parse 3D datagrids
            while endblock not in data[i + j]:
                if datatype in data[i + j] and datatype == "DATAGRID_3D":
                    nvol = list(map(int, data[i + j + 1].split()))
                    vol = [[[0] * nvol[2] for y in range(nvol[1])]
                           for z in range(nvol[0])]
                    origin = list(map(float, data[i + j + 2].split()))
                    # assume spanning space equals cell
                    j += 6
                    line = data[i + j].split()
                    for z in range(nvol[2]):
                        for y in range(nvol[1]):
                            for x in range(nvol[0]):
                                while not (line or endgrid in data[i + j + 1]):
                                    j += 1
                                    line = data[i + j].split()
                                vol[x][y][z] = float(line.pop(0))
                    del vol[-1]
                    for n in vol:
                        del n[-1]
                        for m in n:
                            del m[-1]
                    grids.append((vol, origin))
                    j += 1
                j += 1
            if grids:
                tmol.changeStep(0)
                tmol.setVol(*grids[0])
                if len(grids) > 1:
                    for k in range(1, len(grids)):
                        tmol.copyStep()
                        tmol.setVol(*grids[k])
            i += j
        i += 1
    return tmol, None
