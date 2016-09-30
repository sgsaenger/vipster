from io import StringIO
from vipster.ioplugins import empire
from vipster import Molecule
from test_preamble import *


def test_empire_parse():
    name = 'empire_in'
    target = ["3\n",
              "Hamil=PM3 calc=spt Periodic\n",
              "O   0.0 0.0 0.0\n",
              "H   0.9572 0.0 0.0\n",
              "H   -0.23999 0.926627 0.0\n",
              "\n",
              "1.0 0.0 0.0\n",
              "0.0 2.0 0.0\n",
              "0.0 0.0 3.0\n"]
    Mol, _ = empire.parser(name, target)
    assert Mol.name == name
    assert Mol.nat == 3
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['H', 'O']
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ['O', (0, 0, 0), '0.',
                       [False, False, False], False])
    assert atom_equal(Mol.getAtom(2, fmt='bohr'),
                      ['H', (-0.4535, 1.7510, 0)])
    assert len(Mol.getBonds(1.1)[0]) == 2
    assert vec_equal(Mol.getVec(), ((1, 0, 0), (0, 2, 0), (0, 0, 3)))
    assert float_equal(Mol.getCellDim(fmt='angstrom'), 1)


def test_empire_write():
    target = StringIO("3\n"
                      "Hamil=PM3 calc=spt Periodic\n"
                      "O     0.0000  0.0000  0.0000\n"
                      "H     0.9572  0.0000  0.0000\n"
                      "H    -0.2400  0.9266  0.0000\n"
                      "\n"
                      "1.0000 0.0000 0.0000\n"
                      "0.0000 2.0000 0.0000\n"
                      "0.0000 0.0000 3.0000\n")
    f = StringIO()
    Mol = Molecule()
    Mol.newAtom('O', (0, 0, 0), fmt='angstrom')
    Mol.newAtom('H', (0.9572, 0, 0), fmt='angstrom')
    Mol.newAtom('H', (-0.23999, 0.926627, 0), fmt='angstrom')
    Mol.setCellDim(1, fmt='angstrom')
    Mol.setVec(((1, 0, 0), (0, 2, 0), (0, 0, 3)))
    empire.writer(Mol, f, None)
    assert f.getvalue() == target.getvalue()
