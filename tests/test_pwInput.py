from io import StringIO
from vipster.ioplugins import pwInput
from vipster import Molecule
from test_preamble import *


def test_pwi_parse_scf_ibrav0cdm_gamma_alat():
    name = "pwi_scf_ibrav0cdm_gamma_alat"
    target = ["&control\n",
              " calculation='scf'\n",
              "/\n",
              "&system\n",
              " nat=8\n",
              " ntyp=2\n",
              " ibrav=0\n",
              " celldm(1)=9.448631\n",
              " ecutwfc=30.0\n",
              "/\n",
              "&electrons\n",
              " diagonalization='david'\n",
              " conv_thr=1.0e-8\n",
              "/\n",
              "\n",
              "ATOMIC_SPECIES\n",
              "Na  22.99  Na.uspp736.pbe.UPF\n",
              "Cl  35.45  Cl.uspp736.pbe.UPF\n",
              "\n",
              "ATOMIC_POSITIONS alat\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    0.5000  0.0000  0.0000\n",
              "Na    0.5000  0.5000  0.0000\n",
              "Cl    0.0000  0.5000  0.0000\n",
              "Cl    0.0000  0.0000  0.5000\n",
              "Na    0.5000  0.0000  0.5000\n",
              "Cl    0.5000  0.5000  0.5000\n",
              "Na    0.0000  0.5000  0.5000\n",
              "\n",
              "K_POINTS gamma\n",
              "\n",
              "CELL_PARAMETERS\n",
              "1.00000 0.00000 0.00000\n",
              "0.00000 1.00000 0.00000\n",
              "0.00000 0.00000 1.00000\n"]
    tparam = {"type": "pwi", "name": name,
              "&control": {"calculation": "'scf'"},
              "&system": {"ibrav": "0", "ecutwfc": "30.0"},
              "&electrons": {"diagonalization": "'david'",
                             "conv_thr": "1.0e-8"}}
    Mol, p = pwInput.parser(name, target)
    assert p == tparam
    assert Mol.name == name
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    assert Mol.getFmt() == 'alat'
    assert Mol.getKpoints('active') == 'gamma'
    assert float_equal(Mol.getCellDim(fmt='angstrom'), 5)
    assert vec_equal(Mol.getVec(), ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]


def test_pwi_parse_vcrelax_ibrav8cdm_mpg_angstrom():
    name = "pwi_vcrelax_ibrav8cdm_mpg_angstrom"
    target = ["&control\n",
              " calculation='vc-relax'\n",
              "/\n",
              "&system\n",
              " nat=8\n",
              " ntyp=2\n",
              " ibrav=8\n",
              " celldm(1)=9.448631\n",
              " celldm(2)=1.06\n",
              " celldm(3)=1\n",
              " ecutwfc=30.0\n",
              "/\n",
              "&electrons\n",
              " diagonalization='david'\n",
              " conv_thr=1.0e-8\n",
              "/\n",
              "&ions\n",
              " ion_dynamics='bfgs'\n",
              "/\n",
              "&cell\n",
              " cell_dynamics='bfgs'\n",
              " cell_dofree='all'\n",
              "/\n",
              "\n",
              "ATOMIC_SPECIES\n",
              "Na  22.99  Na.uspp736.pbe.UPF\n",
              "Cl  35.45  Cl.uspp736.pbe.UPF\n",
              "\n",
              "ATOMIC_POSITIONS angstrom\n",
              "Na    0.0000  0.0000  0.0000\n",
              "Cl    2.5000  0.0000  0.0000\n",
              "Na    2.5000  2.5000  0.0000\n",
              "Cl    0.0000  2.5000  0.0000\n",
              "Cl    0.0000  0.0000  2.5000\n",
              "Na    2.5000  0.0000  2.5000\n",
              "Cl    2.5000  2.5000  2.5000\n",
              "Na    0.0000  2.5000  2.5000\n",
              "\n",
              "K_POINTS automatic\n",
              "2 0 0 0 0 0\n"]
    tparam = {"type": "pwi", "name": name,
              "&control": {"calculation": "'vc-relax'"},
              "&system": {"ibrav": "0", "ecutwfc": "30.0"},
              "&electrons": {"diagonalization": "'david'",
                             "conv_thr": "1.0e-8"},
              "&ions": {"ion_dynamics": "'bfgs'"},
              "&cell": {"cell_dynamics": "'bfgs'",
                        "cell_dofree": "'all'"}}
    Mol, p = pwInput.parser(name, target)
    assert p == tparam
    assert Mol.name == name
    assert Mol.nat == 8
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['Na', 'Cl']
    assert Mol.getFmt() == 'angstrom'
    assert Mol.getKpoints('active') == 'mpg'
    assert Mol.getKpoints('mpg') == ['2', '0', '0', '0', '0', '0']
    assert float_equal(Mol.getCellDim(fmt='angstrom'), 5)
    assert vec_equal(Mol.getVec(), ((1, 0, 0), (0, 1.06, 0), (0, 0, 1)))
    assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 0, 0, 4, 0, 0, 0]


def test_pwi_parse_relax_ibrav14abc_discrete_crystal_fixpos():
    name = "pwi_relax_ibrav14abc_discrete_crystal_fixpos"
    target = ["&control\n",
              " calculation='relax'\n",
              "/\n",
              "&system\n",
              " nat=4\n",
              " ntyp=1\n",
              " ibrav=14\n",
              " A=2\n",
              " B=1\n",
              " C=1.5\n",
              " cosAB=0.5\n",
              " cosAC=0.4\n",
              " cosBC=0.3\n",
              " ecutwfc=30.0\n",
              "/\n",
              "&electrons\n",
              " diagonalization='david'\n",
              " conv_thr=1.0e-8\n",
              "/\n",
              "&ions\n",
              " ion_dynamics='bfgs'\n",
              "/\n",
              "\n",
              "ATOMIC_SPECIES\n",
              "C  12.0107  C.uspp736.pbe.UPF\n",
              "/\n",
              "\n",
              "ATOMIC_POSITIONS crystal\n",
              "C 0.0 0.0 0.0 0 1 0\n",
              "C 0.25 0.75 0.5 1 0\n",
              "C 0.5 0.25 0.75\n",
              "C 0.75 0.5 0.25\n",
              "/\n",
              "\n",
              "K_POINTS crystal_b\n",
              "2\n",
              "0.000 0.000 0.000 0.75\n",
              "0.500 0.500 0.500 0.25\n"]
    tparam = {"type": "pwi", "name": name,
              "&control": {"calculation": "'relax'"},
              "&system": {"ibrav": "0", "ecutwfc": "30.0"},
              "&electrons": {"diagonalization": "'david'",
                             "conv_thr": "1.0e-8"},
              "&ions": {"ion_dynamics": "'bfgs'"}}
    Mol, p = pwInput.parser(name, target)
    assert p == tparam
    assert Mol.name == name
    assert Mol.nat == 4
    assert Mol.ntyp == 1
    assert Mol.getTypes() == ['C']
    assert Mol.getFmt() == 'crystal'
    assert Mol.getKpoints('active') == 'discrete'
    assert Mol.getKpoints('discrete') == [['0.000', '0.000', '0.000', '0.75'],
                                          ['0.500', '0.500', '0.500', '0.25']]
    assert float_equal(Mol.getCellDim(fmt='angstrom'), 2)
    assert vec_equal(Mol.getVec(), ((1, 0, 0), (0.5, 0.866025, 0),
                                    (0.6, 0.173205, 1.363818)))
    assert list(map(len, Mol.getBonds(1.1))) == [3, 1, 2, 3, 1, 2, 2, 0]
    assert atom_equal(Mol.getAtom(0, fix=True, fmt='crystal'),
                      ['C', [0, 0, 0], [True, False, True]])
    assert atom_equal(Mol.getAtom(1, fix=True, fmt='crystal'),
                      ['C', [0.25, 0.75, 0.5], [False, True, False]])


def test_pwi_write():
    target = StringIO(
        "&control\n"
        " calculation='scf'\n"
        "/\n"
        "\n"
        "&system\n"
        " nat=2\n"
        " ntyp=1\n"
        " celldm(1)=2.0\n"
        " ibrav=0\n"
        " ecutwfc=30.0\n"
        "/\n"
        "\n"
        "&electrons\n"
        " diagonalization='david'\n"
        " conv_thr=1.0e-8\n"
        "/\n"
        "\n"
        "ATOMIC_SPECIES\n"
        "C    12.0107   C.uspp736.pbe.UPF\n"
        "\n"
        "ATOMIC_POSITIONS crystal\n"
        "C     0.00000  0.00000  0.00000 1 0 1\n"
        "C     0.50000  0.50000  0.50000\n"
        "\n"
        "K_POINTS automatic\n"
        "2    0    0    0    0    0   \n"
        "\n"
        "CELL_PARAMETERS\n"
        " 1.00000  2.00000  0.00000\n"
        " 0.00000  5.00000  6.00000\n"
        " 0.00000  8.00000  9.00000\n")
    Mol = Molecule()
    Mol.setFmt('crystal')
    Mol.setCellDim(2)
    Mol.setVec(((1, 2, 0), (0, 5, 6), (0, 8, 9)))
    Mol.newAtom('C', (0, 0, 0), fix=[False, True])
    Mol.newAtom('C', [0.5, 0.5, 0.5], fmt='crystal')
    Mol.setKpoints('active', 'mpg')
    Mol.setKpoints('mpg', ('2', '0', '0', '0', '0', '0'))
    Param = pwInput.param['default']
    Param['&control']['calculation'] = "'scf'"
    Param['&electrons']['diagonalization'] = "'david'"
    Param['&electrons']['conv_thr'] = "1.0e-8"
    f = StringIO()
    pwInput.writer(Mol, f, Param)
    assert f.getvalue() == target.getvalue()
