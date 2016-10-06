from test_preamble import *
from vipster.ioplugins import lammpsCustom


def test_lammpsCustom_parse_angstrom():
    name = 'lammps_custom_dump_angstrom'
    target = ["ITEM: TIMESTEP\n",
              "0\n",
              "ITEM: NUMBER OF ATOMS\n",
              "8\n",
              "ITEM: BOX BOUNDS ff pp ff\n",
              "0 5\n",
              "0 5\n",
              "0 5\n",
              "ITEM: ATOMS id element x y z\n",
              "1 Na   0.0000  0.0000  0.0000\n",
              "2 Cl   2.5000  0.0000  0.0000\n",
              "3 Na   2.5000  2.5000  0.0000\n",
              "4 Cl   0.0000  2.5000  0.0000\n",
              "5 Cl   0.0000  0.0000  2.5000\n",
              "6 Na   2.5000  0.0000  2.5000\n",
              "7 Cl   2.5000  2.5000  2.5000\n",
              "8 Na   0.0000  2.5000  2.5000\n",
              "\n",
              "ITEM: TIMESTEP\n",
              "1\n",
              "ITEM: NUMBER OF ATOMS\n",
              "8\n",
              "ITEM: BOX BOUNDS ff pp ff\n",
              "0 5\n",
              "0 5.3\n",
              "0 5\n",
              "ITEM: ATOMS id element x y z\n",
              "1 Na   0.0000  0.0000  0.0000\n",
              "2 Cl   2.5000  0.0000  0.0000\n",
              "3 Na   2.5000  2.5000  0.0000\n",
              "4 Cl   0.0000  2.5000  0.0000\n",
              "5 Cl   0.0000  0.0000  2.5000\n",
              "6 Na   2.5000  0.0000  2.5000\n",
              "7 Cl   2.5000  2.5000  2.5000\n",
              "8 Na   0.0000  2.5000  2.5000\n"]
    Mol, _ = lammpsCustom.parser(name, target)
    assert len(Mol) == 2
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    for step in Mol.steps:
        assert step.getFmt() == 'angstrom'
        assert float_equal(step.getCellDim(fmt='angstrom'), 1)
    assert vec_equal(Mol.steps[0].getVec(),
                     ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    assert list(map(len, Mol.steps[0].getBonds(1.1))) ==\
        [12, 4, 4, 0, 4, 0, 0, 0]
    assert vec_equal(Mol.steps[1].getVec(),
                     ((5, 0, 0), (0, 5.3, 0), (0, 0, 5)))
    assert list(map(len, Mol.steps[1].getBonds(1.1))) ==\
        [12, 4, 0, 0, 4, 0, 0, 0]
    for i in range(8):
        assert atom_equal(Mol.steps[0].getAtom(i, fmt='bohr'),
                          Mol.steps[1].getAtom(i, fmt='bohr'))


def test_lammpsCustom_parse_crystal_charge():
    name = 'lammps_custom_dump_crystal_charged'
    target = ["ITEM: TIMESTEP\n",
              "0\n",
              "ITEM: NUMBER OF ATOMS\n",
              "8\n",
              "ITEM: BOX BOUNDS ff pp ff\n",
              "2 7\n",
              "2 7\n",
              "2 7\n",
              "ITEM: ATOMS id element xs ys zs q\n",
              "1 Na   0.0000  0.0000  0.0000  0.5\n",
              "2 Cl   0.5000  0.0000  0.0000 -0.5\n",
              "3 Na   0.5000  0.5000  0.0000  0.5\n",
              "4 Cl   0.0000  0.5000  0.0000 -0.5\n",
              "5 Cl   0.0000  0.0000  0.5000 -0.5\n",
              "6 Na   0.5000  0.0000  0.5000  0.5\n",
              "7 Cl   0.5000  0.5000  0.5000 -0.5\n",
              "8 Na   0.0000  0.5000  0.5000  0.5\n"]
    Mol, _ = lammpsCustom.parser(name, target)
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    assert Mol.getFmt() == 'crystal'
    assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]
    assert float_equal(Mol.getCellDim(fmt='angstrom'), 1)
    assert vec_equal(Mol.getVec(), ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    assert atom_equal(Mol.getAtom(1, charge=True, fmt='angstrom'),
                      ['Cl', [2.5, 0, 0], -0.5])
    assert atom_equal(Mol.getAtom(2, charge=True, fmt='angstrom'),
                      ['Na', [2.5, 2.5, 0], 0.5])


def test_lammpsCustom_parse_invalid():
    name = 'lammps_invalid_dump'
    target = ["ITEM: TIMESTEP\n",
              "0\n",
              "ITEM: NUMBER OF ATOMS\n",
              "8\n",
              "ITEM: BOX BOUNDS ff pp ff\n",
              "2 7\n",
              "2 7\n",
              "2 7\n",
              "ITEM: ATOMS id xs ys zs\n",
              "1 0.0000  0.0000  0.0000\n",
              "2 0.5000  0.0000  0.0000\n",
              "3 0.5000  0.5000  0.0000\n",
              "4 0.0000  0.5000  0.0000\n",
              "5 0.0000  0.0000  0.5000\n",
              "6 0.5000  0.0000  0.5000\n",
              "7 0.5000  0.5000  0.5000\n",
              "8 0.0000  0.5000  0.5000\n"]

    def catchNotImplementedError():
        try:
            lammpsCustom.parser(name, target)
        except NotImplementedError:
            return True
        except:
            return False
    assert catchNotImplementedError()
