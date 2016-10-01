from vipster.ioplugins import aimall
from test_preamble import atom_equal


def test_aimall_parse():
    name = 'aimall'
    target = ["Header\n",
              "Header\n",
              "\n",
              "Header\n",
              "CP# 1      Coords =  0.00000000000000E+00  "
              "0.00000000000000E+00  0.00000000000000E+00\n",
              "           Type = (3,-3) NACP O1\n",
              "           Additional Properties\n",
              "CP# 2      Coords =  1.80890000000000E+00  "
              "0.00000000000000E+00  0.00000000000000E+00\n",
              "           Type = (3,-3) NACP H2\n",
              "           Additional Properties\n",
              "CP# 3      Coords = -0.45350000000000E+00  "
              "1.75110000000000E+00  0.00000000000000E+00\n",
              "           Type = (3,-3) NACP H3\n",
              "           Additional Properties\n",
              "CP# 4      Coords =  1.00000000000000E+00  "
              "0.00000000000000E+00  0.00000000000000E+00\n",
              "           Type = (3,-1) BCP O1 H2\n",
              "           Additional Properties\n"]
    Mol, _ = aimall.parser(name, target)
    assert Mol.name == name
    assert Mol.nat == 4
    assert Mol.ntyp == 3
    assert Mol.getTypes() == ['BCP', 'H', 'O']
    assert Mol.getFmt() == 'bohr'
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ['O', (0, 0, 0), '0.',
                       [False, False, False], False])
    assert atom_equal(Mol.getAtom(2, fmt='bohr'),
                      ['H', (-0.4535, 1.7511, 0)])
    assert atom_equal(Mol.getAtom(1, fmt='angstrom'),
                      ['H', (0.9572, 0, 0)])
    assert len(Mol.getBonds(1.1)[0]) == 2
