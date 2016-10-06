from vipster.ioplugins import pwOutput
from test_preamble import *
# flake8: noqa


def test_pwOut_parse_fixcell_gamma_crystal_fixpos():
    name = 'pwo_fixcell_gamma_crystal_fixpos'
    target = ["Headerstuff\n",
              "   gamma-point specific algorithms are used\n",
              "more stuff\n",
              "   number of atoms/cell    =    8\n",
              "some more stuff\n",
              " celldm(1)=    4.724353  celldm(2)=   0.000000  celldm(3)=  0.00000\n",
              " celldm(4)=    0.000000  celldm(5)=   0.000000  celldm(6)=  0.00000\n",
              "\n",
              " crystal axes: (cart. coord. in units of alat)\n",
              "           a(1) = (   2.000000   0.000000   0.000000\n",
              "           a(2) = (   0.000000   2.000000   0.000000\n",
              "           a(3) = (   0.000000   0.000000   2.000000\n",
              "even more stuff\n",
              "  site n.     atom                   positions (alat units)\n",
              "      1           Na   tau(   1) = (   0.0000000   0.0000000   0.0000000  )\n",
              "      2           Cl   tau(   2) = (   1.0000000   0.0000000   0.0000000  )\n",
              "      3           Na   tau(   3) = (   1.0000000   1.0000000   0.0000000  )\n",
              "      4           Cl   tau(   4) = (   0.0000000   1.0000000   0.0000000  )\n",
              "      5           Cl   tau(   5) = (   0.0000000   0.0000000   1.0000000  )\n",
              "      6           Na   tau(   6) = (   1.0000000   0.0000000   1.0000000  )\n",
              "      7           Cl   tau(   7) = (   1.0000000   1.0000000   1.0000000  )\n",
              "      8           Na   tau(   8) = (   0.0000000   1.0000000   1.0000000  )\n",
              "almost done\n",
              "ATOMIC_POSITIONS (crystal)\n",
              "Na   0.0000  0.0000  0.0000 0 0 0\n",
              "Cl   0.5000  0.0000  0.0000 0\n",
              "Na   0.5000  0.5000  0.0000\n",
              "Cl   0.0000  0.5000  0.0000\n",
              "Cl   0.0000  0.0000  0.5000\n",
              "Na   0.5000  0.0000  0.5000\n",
              "Cl   0.5000  0.5000  0.5000\n",
              "Na   0.0000  0.5000  0.5000\n",
              "last filler\n",
              "Begin final coordinates\n"]
    Mol, _ = pwOutput.parser(name, target)
    assert len(Mol) == 2
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    assert Mol.getKpoints('active') == 'gamma'
    assert Mol.getAtom(0, fix=True)[2] == [True, True, True]
    assert Mol.getAtom(1, fix=True)[2] == [True, False, False]
    assert Mol.getAtom(2, fix=True)[2] == [False, False, False]
    assert Mol.steps[0].getFmt() == 'alat'
    assert Mol.steps[1].getFmt() == 'crystal'
    for step in Mol.steps:
        assert list(map(len, step.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]
        assert float_equal(step.getCellDim(fmt='angstrom'), 2.5)
        assert vec_equal(step.getVec(), ((2, 0, 0), (0, 2, 0), (0, 0, 2)))
    for i in range(8):
        assert atom_equal(Mol.steps[0].getAtom(i, fmt='bohr'),
                          Mol.steps[1].getAtom(i, fmt='bohr'))


def test_pwOut_parse_varcell_discrete_angstrom():
    name = 'pwo_fixcell_gamma_crystal_fixpos'
    target = ["Headerstuff\n",
              "   number of atoms/cell    =    8\n",
              "more stuff\n",
              " celldm(1)=    9.448631  celldm(2)=   0.000000  celldm(3)=  0.00000\n",
              " celldm(4)=    0.000000  celldm(5)=   0.000000  celldm(6)=  0.00000\n",
              "\n",
              " crystal axes: (cart. coord. in units of alat)\n",
              "           a(1) = (   1.000000   0.000000   0.000000\n",
              "           a(2) = (   0.000000   1.000000   0.000000\n",
              "           a(3) = (   0.000000   0.000000   1.000000\n",
              "some more stuff\n",
              "  site n.     atom                   positions (alat units)\n",
              "      1           Na   tau(   1) = (   0.0000000   0.0000000   0.0000000  )\n",
              "      2           Cl   tau(   2) = (   0.5000000   0.0000000   0.0000000  )\n",
              "      3           Na   tau(   3) = (   0.5000000   0.5000000   0.0000000  )\n",
              "      4           Cl   tau(   4) = (   0.0000000   0.5000000   0.0000000  )\n",
              "      5           Cl   tau(   5) = (   0.0000000   0.0000000   0.5000000  )\n",
              "      6           Na   tau(   6) = (   0.5000000   0.0000000   0.5000000  )\n",
              "      7           Cl   tau(   7) = (   0.5000000   0.5000000   0.5000000  )\n",
              "      8           Na   tau(   8) = (   0.0000000   0.5000000   0.5000000  )\n",
              "\n",
              "number of k points=    2 gaussian smearing, width (Ry)=  0.0100\n",
              "                  cart. coord. in units 2pi/alat\n",
              "   k(    1) = (   0.0000000   0.0000000   0.0000000), wk =   0.7500000\n",
              "   k(    2) = (   0.5000000   0.0000000   0.0000000), wk =   0.2500000\n",
              "even more stuff\n",
              "CELL_PARAMETERS (alat=  9.448631)\n",
              "   1.000000000   0.000000000   0.000000000\n",
              "   0.000000000   1.060000000   0.000000000\n",
              "   0.000000000   0.000000000   1.000000000\n",
              "\n",
              "ATOMIC_POSITIONS (angstrom)\n",
              "Na   0.0000  0.0000  0.0000\n",
              "Cl   2.5000  0.0000  0.0000\n",
              "Na   2.5000  2.5000  0.0000\n",
              "Cl   0.0000  2.5000  0.0000\n",
              "Cl   0.0000  0.0000  2.5000\n",
              "Na   2.5000  0.0000  2.5000\n",
              "Cl   2.5000  2.5000  2.5000\n",
              "Na   0.0000  2.5000  2.5000\n",
              "last filler\n",
              "Begin final coordinates\n"]
    Mol, _ = pwOutput.parser(name, target)
    assert len(Mol) == 2
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    assert Mol.getKpoints('active') == 'discrete'
    assert Mol.getKpoints('discrete') ==\
        [['0.0000000', '0.0000000', '0.0000000', '0.7500000'],
         ['0.5000000', '0.0000000', '0.0000000', '0.2500000']]
    assert Mol.steps[0].getFmt() == 'alat'
    assert Mol.steps[1].getFmt() == 'angstrom'
    assert list(map(len, Mol.steps[0].getBonds(1.1))) ==\
        [12, 4, 4, 0, 4, 0, 0, 0]
    assert float_equal(Mol.steps[0].getCellDim(fmt='angstrom'), 5)
    assert vec_equal(Mol.steps[0].getVec(), ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    assert list(map(len, Mol.steps[1].getBonds(1.1))) ==\
        [12, 4, 0, 0, 4, 0, 0, 0]
    assert float_equal(Mol.steps[1].getCellDim(fmt='angstrom'), 5)
    assert vec_equal(Mol.steps[1].getVec(),
                     ((1, 0, 0), (0, 1.06, 0), (0, 0, 1)))
    for i in range(8):
        assert atom_equal(Mol.steps[0].getAtom(i, fmt='bohr'),
                          Mol.steps[1].getAtom(i, fmt='bohr'))
