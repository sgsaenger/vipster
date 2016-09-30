from vipster.ioplugins import xsf
from test_preamble import *


def test_xsf_parse_molecular():
    name = 'xsf_molecular'
    target = ["# comment\n",
              "ATOMS\n",
              "O    0.0000  0.0000  0.0000\n",
              "H    0.9572  0.0000  0.0000\n",
              "H   -0.2400  0.9266  0.0000"]
    Mol, _ = xsf.parser(name, target)
    assert Mol.name == name
    assert Mol.nat == 3
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['H', 'O']
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ['O', (0, 0, 0), '0.',
                       [False, False, False], False])
    assert atom_equal(Mol.getAtom(2, fmt='bohr'),
                      ['H', (-0.4535, 1.7511, 0)])
    assert atom_equal(Mol.getAtom(1, fmt='angstrom'),
                      ['H', (0.9572, 0, 0)])
    assert len(Mol.getBonds(1.1)[0]) == 2


def test_xsf_parse_crystal():
    name = 'xsf_crystal'
    target = ["CRYSTAL\n",
              "PRIMVEC\n",
              "   5.0000  0.0000  0.0000\n",
              "   0.0000  5.0000  0.0000\n",
              "   0.0000  0.0000  5.0000\n",
              "CONVCOORD\n",
              "   5.0000  0.0000  0.0000\n",
              "   0.0000  5.0000  0.0000\n",
              "   0.0000  0.0000  5.0000\n",
              "PRIMCOORD\n",
              "8 1\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.5000  0.0000\n",
              "Cl    0.0000  2.5000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.5000  2.5000\n",
              "Na    0.0000  2.5000  2.5000"]
    Mol, _ = xsf.parser(name, target)
    assert Mol.name == name
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]


def test_xsf_parse_molecule_animated():
    name = 'xsf_molecular_anim'
    target = ["ANIMSTEPS 3\n",
              "ATOMS 1\n",
              "O    0.0000  0.0000  0.0000\n",
              "H    0.9572  0.0000  0.0000\n",
              "H   -0.2400  0.9266  0.0000\n",
              "ATOMS 2\n",
              "O    0.0000  0.0000  0.0000\n",
              "H   -0.9572  0.0000  0.0000\n",
              "H    0.2400  0.9266  0.0000\n",
              "ATOMS 3\n",
              "O    0.0000  0.0000  0.0000\n",
              "H   -0.9572  0.0000  0.0000\n",
              "H    0.2400 -0.9266  0.0000"]
    Mol, _ = xsf.parser(name, target)
    assert len(Mol) == 3
    assert Mol.name == name
    for step in Mol.steps:
        assert step.nat == 3
        assert step.ntyp == 2
        assert step.getTypes() == ['H', 'O']
        assert len(step.getBonds(1.1)[0]) == 2


def test_xsf_parse_crystal_animated_fixcell():
    name = 'xsf_crystal_anim_fixcell'
    target = ["ANIMSTEPS 2\n",
              "CRYSTAL\n",
              "PRIMVEC\n",
              "   5.0000  0.0000  0.0000\n",
              "   0.0000  5.3000  0.0000\n",
              "   0.0000  0.0000  5.0000\n",
              "PRIMCOORD 1\n",
              "8 1\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.5000  0.0000\n",
              "Cl    0.0000  2.5000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.5000  2.5000\n",
              "Na    0.0000  2.5000  2.5000\n",
              "PRIMCOORD 2",
              "8 1\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.6000  0.0000\n",
              "Cl    0.0000  2.6000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.6000  2.5000\n",
              "Na    0.0000  2.6000  2.5000"]
    Mol, _ = xsf.parser(name, target)
    assert Mol.name == name
    for step in Mol.steps:
        assert step.nat == 8
        assert step.ntyp == 2
        assert step.getTypes() == ['Na', 'Cl']
    assert list(map(len, Mol.steps[0].getBonds(1.1))) ==\
        [12, 4, 0, 0, 4, 0, 0, 0]
    assert list(map(len, Mol.steps[1].getBonds(1.1))) ==\
        [12, 4, 4, 0, 4, 0, 0, 0]


def test_xsf_parse_crystal_animated_varcell():
    name = 'xsf_crystal_anim_varcell'
    target = ["ANIMSTEPS 2\n",
              "CRYSTAL\n",
              "PRIMVEC 1\n",
              "   5.0000  0.0000  0.0000\n",
              "   0.0000  5.0000  0.0000\n",
              "   0.0000  0.0000  5.0000\n",
              "PRIMCOORD 1\n",
              "8 1\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.5000  0.0000\n",
              "Cl    0.0000  2.5000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.5000  2.5000\n",
              "Na    0.0000  2.5000  2.5000\n",
              "PRIMVEC 2\n",
              "   5.0000  0.0000  0.0000\n",
              "   0.0000  5.3000  0.0000\n",
              "   0.0000  0.0000  5.0000\n",
              "PRIMCOORD 2",
              "8 1\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.5000  0.0000\n",
              "Cl    0.0000  2.5000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.5000  2.5000\n",
              "Na    0.0000  2.5000  2.5000"]
    Mol, _ = xsf.parser(name, target)
    assert Mol.name == name
    for step in Mol.steps:
        assert step.nat == 8
        assert step.ntyp == 2
        assert step.getTypes() == ['Na', 'Cl']
    assert list(map(len, Mol.steps[0].getBonds(1.1))) ==\
        [12, 4, 4, 0, 4, 0, 0, 0]
    assert list(map(len, Mol.steps[1].getBonds(1.1))) ==\
        [12, 4, 0, 0, 4, 0, 0, 0]


def test_xsf_parse_datagrid3d():
    name = 'xsf_datagrid3d'
    target = ["CRYSTAL\n",
              "PRIMVEC\n",
              "   5.0000  0.0000  0.0000\n",
              "   0.0000  5.0000  0.0000\n",
              "   0.0000  0.0000  5.0000\n",
              "PRIMCOORD\n",
              "8 1\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.5000  0.0000\n",
              "Cl    0.0000  2.5000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.5000  2.5000\n",
              "Na    0.0000  2.5000  2.5000\n",
              "BEGIN_BLOCK_DATAGRID_3D\n",
              "3D_TESTGRID\n",
              "BEGIN_DATAGRID_3D_UNKNOWN1\n",
              "  5  5  5\n",
              "  0.0 0.0 0.0\n",
              "  5.0 0.0 0.0\n",
              "  0.0 5.0 0.0\n",
              "  0.0 0.0 5.0\n",
              "    0.000  1.000  2.000  5.196  8.000\n",
              "    1.000  1.414  2.236  5.292  8.062\n",
              "    2.000  2.236  2.828  5.568  8.246\n",
              "    3.000  3.162  3.606  6.000  8.544\n",
              "    4.000  4.123  4.472  6.557  8.944\n",
              "\n",
              "    1.000  1.414  2.236  5.292  8.062\n",
              "    1.414  1.732  2.449  5.385  8.124\n",
              "    2.236  2.449  3.000  5.657  8.307\n",
              "    3.162  3.317  3.742  6.083  8.602\n",
              "    4.123  4.243  4.583  6.633  9.000\n",
              "\n",
              "    2.000  2.236  2.828  5.568  8.246\n",
              "    2.236  2.449  3.000  5.657  8.307\n",
              "    2.828  3.000  3.464  5.916  8.485\n",
              "    3.606  3.742  4.123  6.325  8.775\n",
              "    4.472  4.583  4.899  6.856  9.165\n",
              "\n",
              "    3.000  3.162  3.606  6.000  8.544\n",
              "    3.162  3.317  3.742  6.083  8.602\n",
              "    3.606  3.742  4.123  6.325  8.775\n",
              "    4.243  4.359  4.690  6.708  9.055\n",
              "    5.000  5.099  5.385  7.211  9.434\n",
              "\n",
              "    4.000  4.123  4.472  6.557  8.944\n",
              "    4.123  4.243  4.583  6.633  9.000\n",
              "    4.472  4.583  4.899  6.856  9.165\n",
              "    5.000  5.099  5.385  7.211  9.434\n",
              "    5.657  5.745  6.000  7.681  9.798\n",
              "END_DATAGRID_3D\n",
              "BEGIN_DATAGRID_3D_UNKNOWN2\n",
              "  3  3  3\n",
              "  0.5 0.5 0.5\n",
              "  5.0 0.0 0.0\n",
              "  0.0 5.0 0.0\n",
              "  0.0 0.0 5.0\n",
              "    0.000  1.000 -3.000  2.000  3.000\n",
              "   -3.000 -2.000 -2.000 -3.000  4.000\n",
              "    5.000 -3.000  6.000  7.000 -3.000\n",
              "   -2.000 -2.000 -3.000 -1.000 -1.000\n",
              "   -1.000 -1.000 -1.000 -1.000 -1.000\n",
              "   -1.000 -1.000\n",
              "END_DATAGRID_3D\n",
              "END_BLOCK_DATAGRID_3D\n"]
    Mol, _ = xsf.parser(name, target)
    assert len(Mol) == 2
    assert vec_equal(Mol.steps[0].getVolOffset(), [0, 0, 0])
    assert vec_equal(Mol.steps[1].getVolOffset(), [0.5, 0.5, 0.5])
    assert vec_equal(Mol.getVol(), (((0, 4), (2, 6)),
                                    ((1, 5), (3, 7))))
    assert vec_equal(Mol.getVolGradient(), 0)
    assert Mol.getVolGradient().shape == (3, 2, 2, 2)
