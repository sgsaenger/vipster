from vipster.ioplugins import mol2
from test_preamble import atom_equal


def test_mol2_parse():
    name = 'mol2'
    target = ["@<TRIPOS>MOLECULE\n",
              "MOL\n",
              "   3   2   1   0   0\n",
              "SMALL\n",
              "resp\n",
              "\n",
              "\n",
              "@<TRIPOS>ATOM\n",
              "    1 O1   0.0000   0.0000   0.0000 o     1 MOL    -0.830000\n",
              "    2 H1   0.9572   0.0000   0.0000 h     1 MOL     0.415000\n",
              "    3 H2  -0.2400   0.9266   0.0000 h     1 MOL     0.415000\n",
              "@<TRIPOS>BOND\n",
              "   1   1   2 1\n",
              "   2   1   3 1\n",
              "@<TRIPOS>SUBSTRUCTURE\n",
              "   1 MOL       1 TEMP      0 ****  ****    0 ROOT"]
    Mol, _ = mol2.parser(name, target)
    print(Mol.getAtom(2, fmt='bohr'))
    assert Mol.name == name
    assert Mol.nat == 3
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['H', 'O']
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ['O', (0, 0, 0), '-0.830000',
                       [False, False, False], False])
    assert atom_equal(Mol.getAtom(2, fmt='bohr'),
                      ['H', (-0.4535, 1.7510, 0)])
    assert len(Mol.getBonds(1.1)[0]) == 2
