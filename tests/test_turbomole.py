from vipster.ioplugins import turbomole
from test_preamble import atom_equal


def test_turbomole_parse():
    name = 'turbomole'
    target = ["$coord\n",
              " 0.0000 0.0000 0.0000 o\n",
              " 1.8089 0.0000 0.0000 h\n",
              "-0.4535 1.7511 0.0000 h\n",
              "$end"]
    Mol, _ = turbomole.parser(name, target)
    assert Mol.name == name
    assert Mol.nat == 3
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['H', 'O']
    assert Mol.getFmt() == 'bohr'
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ['O', (0, 0, 0), '0.',
                       [False, False, False], False])
    assert atom_equal(Mol.getAtom(2, fmt='bohr'),
                      ['H', (-0.4535, 1.7511, 0)])
    assert atom_equal(Mol.getAtom(1, fmt='angstrom'),
                      ['H', (0.9572, 0, 0)])
    assert len(Mol.getBonds(1.1)[0]) == 2
