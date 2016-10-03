from io import StringIO
from vipster.ioplugins import pwInput
from vipster import Molecule, newParam
from test_preamble import *


def test_pwi_parse_formats():
    name = "pwi_parse_formats"
    header = ["&control\n",
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
              "\n"]
    formats = {"bohr": ["ATOMIC_POSITIONS bohr\n",
                        "Na    0.0000  0.0000  0.0000\n",
                        "Cl    4.7243  0.0000  0.0000\n",
                        "Na    4.7243  4.7243  0.0000\n",
                        "\n",
                        "Cl    0.0000  4.7243  0.0000\n",
                        "Cl    0.0000  0.0000  4.7243\n",
                        "Na    4.7243  0.0000  4.7243\n",
                        "Cl    4.7243  4.7243  4.7243\n",
                        "Na    0.0000  4.7243  4.7243\n",
                        "\n"],
               "alat": ["ATOMIC_POSITIONS alat\n",
                        "Na    0.0000  0.0000  0.0000\n",
                        "Cl    0.5000  0.0000  0.0000\n",
                        "Na    0.5000  0.5000  0.0000\n",
                        "Cl    0.0000  0.5000  0.0000\n",
                        "Cl    0.0000  0.0000  0.5000\n",
                        "Na    0.5000  0.0000  0.5000\n",
                        "Cl    0.5000  0.5000  0.5000\n",
                        "Na    0.0000  0.5000  0.5000\n",
                        "\n"],
               "angstrom": ["ATOMIC_POSITIONS angstrom\n",
                            "Na    0.0000  0.0000  0.0000\n",
                            "Cl    2.5000  0.0000  0.0000\n",
                            "Na    2.5000  2.5000  0.0000\n",
                            "Cl    0.0000  2.5000  0.0000\n",
                            "Cl    0.0000  0.0000  2.5000\n",
                            "Na    2.5000  0.0000  2.5000\n",
                            "Cl    2.5000  2.5000  2.5000\n",
                            "Na    0.0000  2.5000  2.5000\n",
                            "\n"],
               "crystal": ["ATOMIC_POSITIONS crystal\n",
                           "Na    0.0000  0.0000  0.0000\n",
                           "Cl    0.5000  0.0000  0.0000\n",
                           "Na    0.5000  0.5000  0.0000\n",
                           "Cl    0.0000  0.5000  0.0000\n",
                           "Cl    0.0000  0.0000  0.5000\n",
                           "Na    0.5000  0.0000  0.5000\n",
                           "Cl    0.5000  0.5000  0.5000\n",
                           "Na    0.0000  0.5000  0.5000\n",
                           "\n"]}
    tail = ["K_POINTS gamma\n",
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
    for f in formats:
        target = header + formats[f] + tail
        Mol, p = pwInput.parser(name, target)
        assert p == tparam
        assert len(Mol) == 1
        assert Mol.name == name
        assert Mol.nat == 8
        assert Mol.ntyp == 2
        assert Mol.getTypes() == ['Na', 'Cl']
        assert Mol.getFmt() == f
        assert Mol.getKpoints('active') == 'gamma'
        assert float_equal(Mol.getCellDim(fmt='angstrom'), 5)
        assert vec_equal(Mol.getVec(), ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
        assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]
        assert all(map(lambda x: atom_equal(*x), zip(
            Mol.getAtoms(fmt='angstrom'),
            (('Na', (0, 0, 0)), ('Cl', (2.5, 0, 0)), ('Na', (2.5, 2.5, 0)),
             ('Cl', (0, 2.5, 0)), ('Cl', (0, 0, 2.5)), ('Na', (2.5, 0, 2.5)),
             ('Cl', (2.5, 2.5, 2.5)), ('Na', (0, 2.5, 2.5))))))


def test_pwi_parse_kpoints():
    name = 'pwi_parse_kpoints'
    body = ["&control\n",
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
            "CELL_PARAMETERS\n",
            "1.00000 0.00000 0.00000\n",
            "0.00000 1.00000 0.00000\n",
            "0.00000 0.00000 1.00000\n",
            "\n"]
    kpoints = [["K_POINTS gamma\n"],
               ["K_POINTS automatic\n",
                "2   0   0   0   0   0   \n"],
               ["K_POINTS tpiba\n",
                "2\n",
                "0.0 0.0 0.0 0.75\n",
                "0.5 0.5 0.5 0.25\n"],
               ["K_POINTS tpiba_b\n",
                "2\n",
                "0.0 0.0 0.0 0.75\n",
                "0.5 0.5 0.5 0.25\n"],
               ["K_POINTS crystal\n",
                "2\n",
                "0.0 0.0 0.0 0.75\n",
                "0.5 0.5 0.5 0.25\n"],
               ["K_POINTS crystal_b\n",
                "2\n",
                "0.0 0.0 0.0 0.75\n",
                "0.5 0.5 0.5 0.25\n"]]
    activetarget = ['gamma', 'mpg', 'discrete',
                    'discrete', 'discrete', 'discrete']
    opttarget = [{'crystal': False, 'bands': False},
                 {'crystal': False, 'bands': True},
                 {'crystal': True, 'bands': False},
                 {'crystal': True, 'bands': True}]
    for k in range(6):
        Mol, p = pwInput.parser(name, body + kpoints[k])
        assert Mol.getKpoints('active') == activetarget[k]
        if Mol.getKpoints('active') == 'mpg':
            assert Mol.getKpoints('mpg') == ['2', '0', '0', '0', '0', '0']
        elif Mol.getKpoints('active') == 'discrete':
            assert Mol.getKpoints('discrete') ==\
                [['0.0', '0.0', '0.0', '0.75'],
                 ['0.5', '0.5', '0.5', '0.25']]
            assert Mol.getKpoints('options') == opttarget[k - 2]


def test_pwi_parse_ibrav():
    name = "pwi_parse_ibrav"
    header = ["&control\n",
              " calculation='scf'\n",
              "/\n",
              "&system\n",
              " nat=4\n",
              " ntyp=1\n",
              " ecutwfc=30.0\n"]
    tail = ["/\n",
            "&electrons\n",
            " diagonalization='david'\n",
            " conv_thr=1.0e-8\n",
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
            "K_POINTS gamma\n"]
    fails = [["ibrav=1\n"],
             ["ibrav=4\n", "celldm(1)=9.448631\n"],
             ["ibrav=0\n"]]
    celld = [["ibrav=1\n", "celldm(1)=9.448631\n"],
             ["ibrav=2\n", "celldm(1)=9.448631\n"],
             ["ibrav=3\n", "celldm(1)=9.448631\n"],
             ["ibrav=4\n", "celldm(1)=9.448631\n", "celldm(3)=1.5\n"],
             ["ibrav=5\n", "celldm(1)=9.448631\n", "celldm(4)=0.5\n"],
             ["ibrav=-5\n", "celldm(1)=9.448631\n", "celldm(4)=0.5\n"],
             ["ibrav=6\n", "celldm(1)=9.448631\n", "celldm(3)=1.5\n"],
             ["ibrav=7\n", "celldm(1)=9.448631\n", "celldm(3)=1.5\n"],
             ["ibrav=8\n", "celldm(1)=9.448631\n",
              "celldm(2)=1.25\n", "celldm(3)=1.5\n"],
             ["ibrav=9\n", "celldm(1)=9.448631\n",
              "celldm(2)=1.25\n", "celldm(3)=1.5\n"],
             ["ibrav=10\n", "celldm(1)=9.448631\n",
              "celldm(2)=1.25\n", "celldm(3)=1.5\n"],
             ["ibrav=11\n", "celldm(1)=9.448631\n",
              "celldm(2)=1.25\n", "celldm(3)=1.5\n"],
             ["ibrav=12\n", "celldm(1)=9.448631\n", "celldm(2)=1.25\n",
              "celldm(3)=1.5\n", "celldm(4)=0.5\n"],
             ["ibrav=-12\n", "celldm(1)=9.448631\n", "celldm(2)=1.25\n",
              "celldm(3)=1.5\n", "celldm(5)=0.75\n"],
             ["ibrav=13\n", "celldm(1)=9.448631\n", "celldm(2)=1.25\n",
              "celldm(3)=1.5\n", "celldm(4)=0.5\n"],
             ["ibrav=14\n", "celldm(1)=9.448631\n", "celldm(2)=1.25\n",
              "celldm(3)=1.5\n", "celldm(4)=0.25\n",
              "celldm(5)=0.75\n", "celldm(6)=0.5\n"]]
    cryst = [["ibrav=1\n", "A=5\n"],
             ["ibrav=2\n", "A=5\n"],
             ["ibrav=3\n", "A=5\n"],
             ["ibrav=4\n", "A=5\n", "C=7.5\n"],
             ["ibrav=5\n", "A=5\n", "cosAB=0.5\n"],
             ["ibrav=-5\n", "A=5\n", "cosAB=0.5\n"],
             ["ibrav=6\n", "A=5\n", "C=7.5\n"],
             ["ibrav=7\n", "A=5\n", "C=7.5\n"],
             ["ibrav=8\n", "A=5\n", "B=6.25\n", "C=7.5\n"],
             ["ibrav=9\n", "A=5\n", "B=6.25\n", "C=7.5\n"],
             ["ibrav=10\n", "A=5\n", "B=6.25\n", "C=7.5\n"],
             ["ibrav=11\n", "A=5\n", "B=6.25\n", "C=7.5\n"],
             ["ibrav=12\n", "A=5\n", "B=6.25\n", "C=7.5\n", "cosAB=0.5\n"],
             ["ibrav=-12\n", "A=5\n", "B=6.25\n", "C=7.5\n", "cosAC=0.75\n"],
             ["ibrav=13\n", "A=5\n", "B=6.25\n", "C=7.5\n", "cosAB=0.5\n"],
             ["ibrav=14\n", "A=5\n", "B=6.25\n", "C=7.5\n",
              "cosAB=0.5\n", "cosAC=0.75\n", "cosBC=0.25\n"]]
    vectarget = [((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 ((-0.5, 0, 0.5), (0, 0.5, 0.5), (-0.5, 0.5, 0)),
                 ((0.5, 0.5, 0.5), (-0.5, 0.5, 0.5), (-0.5, -0.5, 0.5)),
                 ((1, 0, 0), (-0.5, np.sqrt(3) * 0.5, 0), (0, 0, 1.5)),
                 ((0.5, -0.288675, 0.81650), (0, 0.57735, 0.81650),
                  (-0.5, -0.288675, 0.81650)),
                 ((0, 0.70711, 0.70711), (0.70711, 0, 0.70711),
                  (0.70711, 0.70711, 0)),
                 ((1, 0, 0), (0, 1, 0), (0, 0, 1.5)),
                 ((0.5, -0.5, 0.75), (0.5, 0.5, 0.75), (-0.5, -0.5, 0.75)),
                 ((1, 0, 0), (0, 1.25, 0), (0, 0, 1.5)),
                 ((0.5, 0.625, 0), (-0.5, 0.625, 0), (0, 0, 1.5)),
                 ((0.5, 0, 0.75), (0.5, 0.625, 0), (0, 0.625, 0.75)),
                 ((0.5, 0.625, 0.75), (-0.5, 0.625, 0.75),
                  (-0.5, -0.625, 0.75)),
                 ((1, 0, 0), (0.625, 1.08253, 0), (0, 0, 1.5)),
                 ((1, 0, 0), (0, 1.25, 0), (1.125, 0, 0.992157)),
                 ((0.5, 0, -0.75), (0.625, 1.08253, 0), (0.5, 0, 0.75)),
                 ((1, 0, 0), (0.625, 1.082532, 0),
                  (1.125, -0.216506, 0.96825))]
    for i in range(len(cryst)):
        MolCell, _ = pwInput.parser(name, header + celld[i] + tail)
        MolCrys, _ = pwInput.parser(name, header + cryst[i] + tail)
        assert float_equal(MolCell.getCellDim(fmt='angstrom'), 5)
        assert vec_equal(MolCell.getVec(), vectarget[i])
        assert float_equal(MolCell.getCellDim(), MolCrys.getCellDim())
        assert vec_equal(MolCell.getVec(), MolCrys.getVec())

    def failPwi(f):
        try:
            MolFail, _ = pwInput.parser(name, header + f + tail)
        except ValueError:
            return True
        else:
            return False
    for f in fails:
        assert failPwi(f)


def test_pwi_write_calc():
    calcs = ["'scf'",
             "'relax'",
             "'vc-relax'"]
    adnl = ["",
            "&ions\n ion_dynamics='bfgs'\n/\n\n",
            "&ions\n ion_dynamics='bfgs'\n/\n\n"
            "&cell\n cell_dynamics='bfgs'\n/\n\n"]
    head = "&control\n calculation="
    body = "\n/\n\n"\
           "&system\n"\
           " nat=2\n"\
           " ntyp=1\n"\
           " celldm(1)=2.0\n"\
           " ibrav=0\n"\
           " ecutwfc=30.0\n"\
           "/\n"\
           "\n"\
           "&electrons\n"\
           " diagonalization='david'\n"\
           " conv_thr=1.0e-8\n"\
           "/\n"\
           "\n"
    tail = "ATOMIC_SPECIES\n"\
           "C    12.0107   C.uspp736.pbe.UPF\n"\
           "\n"\
           "ATOMIC_POSITIONS crystal\n"\
           "C     0.00000  0.00000  0.00000 1 0 1\n"\
           "C     0.50000  0.50000  0.50000\n"\
           "\n"\
           "K_POINTS automatic\n"\
           "2    0    0    0    0    0   \n"\
           "\n"\
           "CELL_PARAMETERS\n"\
           " 1.00000  2.00000  0.00000\n"\
           " 0.00000  5.00000  6.00000\n"\
           " 0.00000  8.00000  9.00000\n"
    Mol = Molecule()
    Mol.setFmt('crystal')
    Mol.setCellDim(2)
    Mol.setVec(((1, 2, 0), (0, 5, 6), (0, 8, 9)))
    Mol.newAtom('C', (0, 0, 0), fix=[False, True])
    Mol.newAtom('C', [0.5, 0.5, 0.5], fmt='crystal')
    Mol.setKpoints('active', 'mpg')
    Mol.setKpoints('mpg', ('2', '0', '0', '0', '0', '0'))

    def failWrite(p):
        try:
            pwInput.writer(Mol, f, p)
        except KeyError:
            return True
        else:
            return False

    for c in range(3):
        for a in range(3):
            Param = newParam('pwi')
            Param['&control']['calculation'] = calcs[c]
            Param['&electrons']['diagonalization'] = "'david'"
            Param['&electrons']['conv_thr'] = "1.0e-8"
            if a > 0:
                Param['&ions'] = {'ion_dynamics': "'bfgs'"}
            if a > 1:
                Param['&cell'] = {'cell_dynamics': "'bfgs'"}
            f = StringIO()
            if a < c:
                failWrite(Param)
            else:
                pwInput.writer(Mol, f, Param)
                target = head + calcs[c] + body + adnl[a] + tail
                assert f.getvalue() == target
    # default to scf when nothing present:
    Param = newParam('pwi')
    Param['&electrons']['diagonalization'] = "'david'"
    Param['&electrons']['conv_thr'] = "1.0e-8"
    f = StringIO()
    pwInput.writer(Mol, f, Param)
    target = "&control" + body + tail
    assert f.getvalue() == target


def test_pwi_write_kpoints():
    body = "&control\n"\
           "/\n"\
           "\n"\
           "&system\n"\
           " nat=2\n"\
           " ntyp=1\n"\
           " celldm(1)=2.0\n"\
           " ibrav=0\n"\
           " ecutwfc=30.0\n"\
           "/\n"\
           "\n"\
           "&electrons\n"\
           "/\n"\
           "\n"\
           "ATOMIC_SPECIES\n"\
           "C    12.0107   C.uspp736.pbe.UPF\n"\
           "\n"\
           "ATOMIC_POSITIONS alat\n"\
           "C     0.00000  0.00000  0.00000 1 0 1\n"\
           "C     0.50000  7.50000  7.50000\n"\
           "\n"
    tail = "\n"\
           "CELL_PARAMETERS\n"\
           " 1.00000  2.00000  0.00000\n"\
           " 0.00000  5.00000  6.00000\n"\
           " 0.00000  8.00000  9.00000\n"
    Param = pwInput.param['default'].copy()
    Mol = Molecule()
    Mol.setFmt('alat')
    Mol.setCellDim(2)
    Mol.setVec(((1, 2, 0), (0, 5, 6), (0, 8, 9)))
    Mol.newAtom('C', (0, 0, 0), fix=[False, True])
    Mol.newAtom('C', [0.5, 0.5, 0.5], fmt='crystal')

    def assertKpoint():
        f = StringIO()
        pwInput.writer(Mol, f, Param)
        for i in range(len(target.split('\n'))):
            assert f.getvalue().split('\n')[i] == target.split('\n')[i]
        assert f.getvalue() == target
    Mol.setKpoints('active', 'gamma')
    target = body + "K_POINTS gamma\n" + tail
    assertKpoint()
    Mol.setKpoints('active', 'mpg')
    Mol.setKpoints('mpg', ('2', '0', '0', '0', '0', '0'))
    target = body +\
        "K_POINTS automatic\n2    0    0    0    0    0   \n" +\
        tail
    assertKpoint()
    Mol.setKpoints('active', 'discrete')
    Mol.setKpoints('discrete',
                   [['0.0', '0.0', '0.0', '0.75'],
                    ['0.5', '0.5', '0.5', '0.25']])
    Mol.setKpoints('options', {'crystal': False, 'bands': False})
    target = body +\
        "K_POINTS tpiba\n2\n0.0  0.0  0.0  0.75\n0.5  0.5  0.5  0.25\n" +\
        tail
    assertKpoint()
    Mol.setKpoints('options', {'crystal': False, 'bands': True})
    target = body +\
        "K_POINTS tpiba_b\n2\n0.0  0.0  0.0  0.75\n0.5  0.5  0.5  0.25\n" +\
        tail
    assertKpoint()
    Mol.setKpoints('options', {'crystal': True, 'bands': False})
    target = body +\
        "K_POINTS crystal\n2\n0.0  0.0  0.0  0.75\n0.5  0.5  0.5  0.25\n" +\
        tail
    assertKpoint()
    Mol.setKpoints('options', {'crystal': True, 'bands': True})
    target = body +\
        "K_POINTS crystal_b\n2\n0.0  0.0  0.0  0.75\n0.5  0.5  0.5  0.25\n" +\
        tail
    assertKpoint()
