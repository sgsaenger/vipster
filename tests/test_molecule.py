from test_preamble import *
from vipster import Molecule


def test_atoms_bohr():
    Mol = Molecule()
    Mol.newAtom()
    Mol.newAtom('H', (1, 1, 1), 1, [True, False], True)
    assert Mol.nat == 2
    assert Mol.ntyp == 2
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ('C', (0, 0, 0), 0., [False, False, False], False))
    assert atom_equal(Mol.getAtom(1, charge=True, fix=True, hidden=True),
                      ('H', (1, 1, 1), 1., [True, False, False], True))
