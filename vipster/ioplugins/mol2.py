# -*- coding: utf-8 -*-
from vipster.molecule import Molecule

name = 'mol2'
extension = 'mol2'
argument = 'mol2'
param = None

writer = None


def parser(name, data):
    """ Parse sybyl mol2 files """
    tmol = Molecule(name)
    parsing = False
    for line in data:
        if not line.strip():
            continue
        elif line.strip()[0] == '#':
            continue
        elif "@<TRIPOS>ATOM" in line:
            parsing = True
        elif parsing:
            if "@<TRIPOS>" in line:
                tmol.setFmt("angstrom", scale=True)
                return tmol, None
            at = line.strip().split()
            tmol.newAtom(at[5].upper(), at[2:5], charge=at[8])
