from vipster.ioplugins import cpmd
from vipster import Molecule, newParam
from io import StringIO
from test_preamble import *


def test_cpmd_parse_format():
    name = 'cpmd_format'
    head = ["&INFO\n",
            "comment\n",
            "&END\n",
            "\n",
            "!comment\n",
            "&CPMD\n",
            " MOLECULAR DYNAMICS CP\n",
            " RESTART WAVEFUNCTION COORDINATES LATEST\n",
            "&END",
            "&SYSTEM\n",
            "CUTOFF\n",
            " 25.0\n",
            "TESR\n",
            " 2\n"]
    tails = [["ANGSTROM\n",
              "CELL VECTORS\n", "5.0000 0.0000 0.0000\n",
              "0.0000 5.0000 0.0000\n", "0.0000 0.0000 5.0000\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 2.50000  2.50000  0.00000\n",
              " 0.00000  2.50000  2.50000\n", " 2.50000  0.00000  2.50000\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  2.50000  0.00000\n", " 2.50000  0.00000  0.00000\n",
              " 0.00000  0.00000  2.50000\n", " 2.50000  2.50000  2.50000\n",
              "&END\n"],
             ["CELL VECTORS\n", "9.4486 0.0000 0.0000\n",
              "0.0000 9.4486 0.0000\n", "0.0000 0.0000 9.4486\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 4.72432  4.72432  0.00000\n",
              " 0.00000  4.72432  4.72432\n", " 4.72432  0.00000  4.72432\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  4.72432  0.00000\n", " 4.72432  0.00000  0.00000\n",
              " 0.00000  0.00000  4.72432\n", " 4.72432  4.72432  4.72432\n",
              "&END\n"],
             ["SCALE CARTESIAN\n",
              "CELL VECTORS\n", "9.4486 0.0000 0.0000\n",
              "0.0000 9.4486 0.0000\n", "0.0000 0.0000 9.4486\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 0.50000  0.50000  0.00000\n",
              " 0.00000  0.50000  0.50000\n", " 0.50000  0.00000  0.50000\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.50000  0.00000\n", " 0.50000  0.00000  0.00000\n",
              " 0.00000  0.00000  0.50000\n", " 0.50000  0.50000  0.50000\n",
              "&END\n"],
             ["SCALE\n",
              "CELL VECTORS\n", "9.4486 0.0000 0.0000\n",
              "0.0000 9.4486 0.0000\n", "0.0000 0.0000 9.4486\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 0.50000  0.50000  0.00000\n",
              " 0.00000  0.50000  0.50000\n", " 0.50000  0.00000  0.50000\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.50000  0.00000\n", " 0.50000  0.00000  0.00000\n",
              " 0.00000  0.00000  0.50000\n", " 0.50000  0.50000  0.50000\n",
              "&END\n"],
             ["SCALE S=0.5\n",
              "CELL VECTORS\n", "9.4486 0.0000 0.0000\n",
              "0.0000 9.4486 0.0000\n", "0.0000 0.0000 9.4486\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 0.25000  0.25000  0.00000\n",
              " 0.00000  0.25000  0.25000\n", " 0.25000  0.00000  0.25000\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.25000  0.00000\n", " 0.25000  0.00000  0.00000\n",
              " 0.00000  0.00000  0.25000\n", " 0.25000  0.25000  0.25000\n",
              "&END\n"],
             ["SCALE SX=0.5\n",
              "CELL VECTORS\n", "9.4486 0.0000 0.0000\n",
              "0.0000 9.4486 0.0000\n", "0.0000 0.0000 9.4486\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 0.25000  0.50000  0.00000\n",
              " 0.00000  0.50000  0.50000\n", " 0.25000  0.00000  0.50000\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.50000  0.00000\n", " 0.25000  0.00000  0.00000\n",
              " 0.00000  0.00000  0.50000\n", " 0.25000  0.50000  0.50000\n",
              "&END\n"],
             ["SCALE SX=0.5 SY=0.5 SZ=2\n",
              "CELL VECTORS\n", "9.4486 0.0000 0.0000\n",
              "0.0000 9.4486 0.0000\n", "0.0000 0.0000 9.4486\n",
              "&END\n",
              "&ATOMS\n",
              "*Na.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.00000  0.00000\n", " 0.25000  0.25000  0.00000\n",
              " 0.00000  0.25000  1.00000\n", " 0.25000  0.00000  1.00000\n",
              "*Cl.uspp.pbe BINARY\n", " LMAX=F\n", " 4\n",
              " 0.00000  0.25000  0.00000\n", " 0.25000  0.00000  0.00000\n",
              " 0.00000  0.00000  1.00000\n", " 0.25000  0.25000  1.00000\n",
              "&END\n"]]
    formats = ['angstrom', 'bohr', 'alat', 'crystal',
               'crystal', 'crystal', 'crystal']
    targetp = newParam('cpmd')
    targetp["name"] = name
    targetp["&INFO"] = "comment\n"
    targetp["&CPMD"] = " MOLECULAR DYNAMICS CP\n RESTART "\
        "WAVEFUNCTION COORDINATES LATEST\n"
    targetp["&SYSTEM"] = "CUTOFF\n 25.0\nTESR\n 2\n"
    for f in range(len(formats)):
        Mol, p = cpmd.parser(name, head + tails[f])
        assert p == targetp
        assert len(Mol) == 1
        assert Mol.name == name
        assert Mol.nat == 8
        assert Mol.ntyp == 2
        assert Mol.getTypes() == ["Na", "Cl"]
        assert Mol.getFmt() == formats[f]
        assert Mol.getKpoints('active') == 'gamma'
        assert vec_equal(Mol.getVec() * Mol.getCellDim(fmt='angstrom'),
                         ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
        assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]
        assert all(map(lambda x: atom_equal(*x), zip(
            Mol.getAtoms(fmt='angstrom'),
            (('Na', (0, 0, 0)), ('Na', (2.5, 2.5, 0)), ('Na', (0, 2.5, 2.5)),
             ('Na', (2.5, 0, 2.5)), ('Cl', (0, 2.5, 0)), ('Cl', (2.5, 0, 0)),
             ('Cl', (0, 0, 2.5)), ('Cl', (2.5, 2.5, 2.5))))))


def test_cpmd_parse_kpoints():
    name = 'cpmd_kpoints'
    head = ["&INFO\n",
            "comment\n",
            "&END\n",
            "\n",
            "!comment\n",
            "&CPMD\n",
            " MOLECULAR DYNAMICS CP\n",
            " RESTART WAVEFUNCTION COORDINATES LATEST\n",
            "&END",
            "&SYSTEM\n"]
    tail = ["CUTOFF\n",
            " 25.0\n",
            "TESR\n",
            " 2\n",
            "SCALE\n",
            "CELL VECTORS\n",
            "9.4486 0.0000 0.0000\n",
            "0.0000 9.4486 0.0000\n",
            "0.0000 0.0000 9.4486\n",
            "&END\n",
            "&ATOMS\n",
            "*Na.uspp.pbe BINARY\n",
            " LMAX=F\n",
            " 4\n",
            " 0.00000  0.00000  0.00000\n",
            " 0.50000  0.50000  0.00000\n",
            " 0.00000  0.50000  0.50000\n",
            " 0.50000  0.00000  0.50000\n",
            "*Cl.uspp.pbe BINARY\n",
            " LMAX=F\n",
            " 4\n",
            " 0.00000  0.50000  0.00000\n",
            " 0.50000  0.00000  0.00000\n",
            " 0.00000  0.00000  0.50000\n",
            " 0.50000  0.50000  0.50000\n",
            "&END\n"]
    ksys = [[],
            ["KPOINTS MONKHORST-PACK\n", "2 0 0\n"],
            ["KPOINTS MONKHORST-PACK SHIFT=0 0 0\n", "2 0 0\n"],
            ["KPOINTS\n", "2\n", "0.0 0.0 0.0 2\n", "0.5 0.5 0.5 0\n"],
            ["KPOINTS SCALED\n", "2\n", "0.0 0.0 0.0 2\n", "0.5 0.5 0.5 0\n"],
            ["KPOINTS BANDS\n", "2 0.0 0.0 0.0 0.5 0.5 0.5\n",
             "0 0 0 0 0 0 0\n"],
            ["KPOINTS SCALED BANDS\n", "2 0.0 0.0 0.0 0.5 0.5 0.5\n",
             "0 0 0 0 0 0 0\n"]]
    kfmt = ['gamma', 'mpg', 'mpg', 'discrete',
            'discrete', 'discrete', 'discrete']
    opttarget = [{'crystal': False, 'bands': False},
                 {'crystal': True, 'bands': False},
                 {'crystal': False, 'bands': True},
                 {'crystal': True, 'bands': True}]
    targetp = newParam('cpmd')
    targetp["name"] = name
    targetp["&INFO"] = "comment\n"
    targetp["&CPMD"] = " MOLECULAR DYNAMICS CP\n RESTART "\
        "WAVEFUNCTION COORDINATES LATEST\n"
    targetp["&SYSTEM"] = "CUTOFF\n 25.0\nTESR\n 2\n"
    for k in range(len(kfmt)):
        Mol, p = cpmd.parser(name, head + ksys[k] + tail)
        assert p == targetp
        assert len(Mol) == 1
        assert Mol.name == name
        assert Mol.nat == 8
        assert Mol.ntyp == 2
        assert Mol.getTypes() == ["Na", "Cl"]
        assert Mol.getFmt() == 'crystal'
        assert vec_equal(Mol.getVec() * Mol.getCellDim(fmt='angstrom'),
                         ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
        assert list(map(len, Mol.getBonds(1.1))) == [12, 4, 4, 0, 4, 0, 0, 0]
        assert all(map(lambda x: atom_equal(*x), zip(
            Mol.getAtoms(fmt='angstrom'),
            (('Na', (0, 0, 0)), ('Na', (2.5, 2.5, 0)), ('Na', (0, 2.5, 2.5)),
             ('Na', (2.5, 0, 2.5)), ('Cl', (0, 2.5, 0)), ('Cl', (2.5, 0, 0)),
             ('Cl', (0, 0, 2.5)), ('Cl', (2.5, 2.5, 2.5))))))
        assert Mol.getKpoints('active') == kfmt[k]
        if Mol.getKpoints('active') == 'mpg':
            assert Mol.getKpoints('mpg') == ['2', '0', '0', '0', '0', '0']
        elif Mol.getKpoints('active') == 'discrete':
            assert Mol.getKpoints('discrete') ==\
                [['0.0', '0.0', '0.0', '2'],
                 ['0.5', '0.5', '0.5', '0']]
            assert Mol.getKpoints('options') == opttarget[k - 3]


def test_cpmd_parse_symmetry():
    name = 'cpmd_symmetry'
    symm = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '12', '14',
            "ISOLATED", "CUBIC", "FACE CENTERED CUBIC", "FCC",
            "BODY CENTERED CUBIC", "BCC", "HEXAGONAL", "TRIGONAL",
            "RHOMBOHEDRAL", "TETRAGONAL", "BODY CENTERED TETRAGONAL",
            "BCT", "ORTHORHOMBIC", "MONOCLINIC", "TRICLINIC", '14', '14', '14']
    cell = [["CELL\n", "5 1 1 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 0 0 0.5 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 1.25 1.5 0 0 0\n"],
            ["CELL\n", "5 1.25 1.5 0.5 0 0\n"],
            ["CELL\n", "5 1.25 1.5 0.25 0.75 0.5\n"],
            ["CELL\n", "5 1 1 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 0 0 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 0 0 0.5 0 0\n"],
            ["CELL\n", "5 0 0 0.5 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 0 1.5 0 0 0\n"],
            ["CELL\n", "5 1.25 1.5 0 0 0\n"],
            ["CELL\n", "5 1.25 1.5 0.5 0 0\n"],
            ["CELL\n", "5 1.25 1.5 0.25 0.75 0.5\n"],
            ["CELL ABSOLUTE\n", "5 6.25 7.5 0.25 0.75 0.5\n"],
            ["CELL DEGREE\n", "5 1.25 1.5 75.5225 41.4096 60\n"],
            ["CELL ABSOLUTE DEGREE\n", "5 6.25 7.5 75.5225 41.4096 60\n"]]
    head = ["&INFO\n",
            "comment\n",
            "&END\n",
            "\n",
            "!comment\n",
            "&CPMD\n",
            " MOLECULAR DYNAMICS CP\n",
            " RESTART WAVEFUNCTION COORDINATES LATEST\n",
            "&END",
            "&SYSTEM\n",
            "CUTOFF\n",
            " 25.0\n",
            "TESR\n",
            " 2\n",
            "ANGSTROM\n",
            "SYMMETRY\n"]
    tail = ["&END\n",
            "&ATOMS\n",
            "*Na.uspp.pbe BINARY\n",
            " LMAX=F\n",
            " 4\n",
            " 0.00000  0.00000  0.00000\n",
            " 0.50000  0.50000  0.00000\n",
            " 0.00000  0.50000  0.50000\n",
            " 0.50000  0.00000  0.50000\n",
            "*Cl.uspp.pbe BINARY\n",
            " LMAX=F\n",
            " 4\n",
            " 0.00000  0.50000  0.00000\n",
            " 0.50000  0.00000  0.00000\n",
            " 0.00000  0.00000  0.50000\n",
            " 0.50000  0.50000  0.50000\n",
            "&END\n"]
    vectarget = [((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 ((-0.5, 0, 0.5), (0, 0.5, 0.5), (-0.5, 0.5, 0)),
                 ((0.5, 0.5, 0.5), (-0.5, 0.5, 0.5), (-0.5, -0.5, 0.5)),
                 ((1, 0, 0), (-0.5, np.sqrt(3) * 0.5, 0), (0, 0, 1.5)),
                 ((0.5, -0.288675, 0.81650), (0, 0.57735, 0.81650),
                  (-0.5, -0.288675, 0.81650)),
                 ((1, 0, 0), (0, 1, 0), (0, 0, 1.5)),
                 ((0.5, -0.5, 0.75), (0.5, 0.5, 0.75), (-0.5, -0.5, 0.75)),
                 ((1, 0, 0), (0, 1.25, 0), (0, 0, 1.5)),
                 ((1, 0, 0), (0.625, 1.08253, 0), (0, 0, 1.5)),
                 ((1, 0, 0), (0.625, 1.082532, 0),
                  (1.125, -0.216506, 0.96825)),
                 ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 ((-0.5, 0, 0.5), (0, 0.5, 0.5), (-0.5, 0.5, 0)),
                 ((-0.5, 0, 0.5), (0, 0.5, 0.5), (-0.5, 0.5, 0)),
                 ((0.5, 0.5, 0.5), (-0.5, 0.5, 0.5), (-0.5, -0.5, 0.5)),
                 ((0.5, 0.5, 0.5), (-0.5, 0.5, 0.5), (-0.5, -0.5, 0.5)),
                 ((1, 0, 0), (-0.5, np.sqrt(3) * 0.5, 0), (0, 0, 1.5)),
                 ((0.5, -0.288675, 0.81650), (0, 0.57735, 0.81650),
                  (-0.5, -0.288675, 0.81650)),
                 ((0.5, -0.288675, 0.81650), (0, 0.57735, 0.81650),
                  (-0.5, -0.288675, 0.81650)),
                 ((1, 0, 0), (0, 1, 0), (0, 0, 1.5)),
                 ((0.5, -0.5, 0.75), (0.5, 0.5, 0.75), (-0.5, -0.5, 0.75)),
                 ((0.5, -0.5, 0.75), (0.5, 0.5, 0.75), (-0.5, -0.5, 0.75)),
                 ((1, 0, 0), (0, 1.25, 0), (0, 0, 1.5)),
                 ((1, 0, 0), (0.625, 1.08253, 0), (0, 0, 1.5)),
                 ((1, 0, 0), (0.625, 1.082532, 0),
                  (1.125, -0.216506, 0.96825)),
                 ((1, 0, 0), (0.625, 1.082532, 0),
                  (1.125, -0.216506, 0.96825)),
                 ((1, 0, 0), (0.625, 1.082532, 0),
                  (1.125, -0.216506, 0.96825)),
                 ((1, 0, 0), (0.625, 1.082532, 0),
                  (1.125, -0.216506, 0.96825))]
    for i in range(len(symm)):
        Mol, _ = cpmd.parser(name, head + [symm[i]] + cell[i] + tail)
        assert float_equal(Mol.getCellDim(fmt='angstrom'), 5)
        assert vec_equal(Mol.getVec(), vectarget[i])


def test_cpmd_parse_constraints():
    name = 'cpmd_constraints'
    head = ["&INFO\n",
            "comment\n",
            "&END\n",
            "\n",
            "!comment\n",
            "&CPMD\n",
            " MOLECULAR DYNAMICS CP\n",
            " RESTART WAVEFUNCTION COORDINATES LATEST\n",
            "&END",
            "&SYSTEM\n",
            "CUTOFF\n",
            " 25.0\n",
            "TESR\n",
            " 2\n",
            "SCALE\n",
            "CELL VECTORS\n",
            "9.4486 0.0000 0.0000\n",
            "0.0000 9.4486 0.0000\n",
            "0.0000 0.0000 9.4486\n",
            "&END\n",
            "&ATOMS\n",
            "*Na.uspp.pbe BINARY\n",
            " LMAX=F\n",
            " 1\n",
            " 0.00000  0.00000  0.00000\n",
            "*Na_MT_PBE.psp KLEINMAN-BYLANDER\n",
            " LMAX=P\n",
            " 1\n",
            " 0.50000  0.50000  0.00000\n",
            "*Na.uspp.pbe BINARY\n",
            " LMAX=F\n",
            " 2\n",
            " 0.00000  0.50000  0.00000\n",
            " 0.50000  0.00000  0.00000\n",
            "ISOTOPE\n",
            " 21.994\n",
            " 22.989\n",
            " 23.990\n",
            "CONSTRAINTS\n"]
    cons = [[],
            ["FIX ALL\n"],
            ["FIX QM\n"],
            ["FIX MM\n"],
            ["FIX PPTY\n", "3\n"],
            ["FIX PPTY SEQUENCE\n", "3 1 3\n"],
            ["FIX ELEMENT\n", "11\n"],
            ["FIX ELEMENT SEQUENCE\n", "11 1 3\n"],
            ["FIX SEQUENCE\n", "1 3\n"],
            ["FIX ATOMS\n", "2\n", "1 3\n"],
            ["FIX COORDINATES\n", "3\n",
             "1 0 1 1\n", "2 1 0 1\n", "3 1 1 0\n"],
            ["FIX COORDINATES\n", "4\n",
             "1 0 0 1\n", "2 0 1 0\n", "3 1 0 0\n", "4 0 0 0\n"],
            ]
    F = [False, False, False]
    X = [True, False, False]
    Y = [False, True, False]
    Z = [False, False, True]
    XY = [True, True, False]
    XZ = [True, False, True]
    YZ = [False, True, True]
    A = [True, True, True]
    fixtarget = [[F, F, F, F], [A, A, A, A], [A, A, A, A], [F, F, F, F],
                 [F, F, A, A], [F, F, A, F], [A, A, A, A], [A, A, A, F],
                 [A, A, A, F], [A, F, A, F], [X, Y, Z, F], [XY, XZ, YZ, A],
                 ]
    for i in range(len(cons)):
        Mol, _ = cpmd.parser(name,
                             head + cons[i] + ["END CONSTRAINTS\n", "&END\n"])
        assert Mol.nat == 4
        assert Mol.ntyp == 3
        assert Mol.getTypes() == ['Na', 'Na1', 'Na2']
        assert Mol.pse['Na']['m'] == 21.994
        assert Mol.pse['Na1']['m'] == 22.989
        assert Mol.pse['Na2']['m'] == 23.990
        assert list(map(lambda x: x[-1], Mol.getAtoms(fix=True))) ==\
            fixtarget[i]


def test_cpmd_write_format():
    Mol = Molecule()
    Mol.setCellDim(1, fmt='angstrom')
    Mol.setVec(((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    Mol.newAtom('Na', (0, 0, 0), fmt='crystal')
    Mol.newAtom('Cl', (0.5, 0.5, 0.5), fmt='crystal')
    p = newParam('cpmd')
    p["&INFO"] = "comment\n"
    p["&CPMD"] = " MOLECULAR DYNAMICS CP\n RESTART "\
        "WAVEFUNCTION COORDINATES LATEST\n"
    fmts = ['bohr', 'angstrom', 'crystal', 'alat']
    head = "&INFO\n"\
           "comment\n"\
           "&END\n"\
           "&CPMD\n"\
           " MOLECULAR DYNAMICS CP\n RESTART "\
           "WAVEFUNCTION COORDINATES LATEST\n"\
           "&END\n"\
           "&SYSTEM\n"
    body = ["  CELL VECTORS\n" "   9.4486  0.0000  0.0000\n"
            "   0.0000  9.4486  0.0000\n" "   0.0000  0.0000  9.4486\n"
            "&END\n" "&ATOMS\n"
            "*Na.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   0.0000  0.0000  0.0000\n"
            "*Cl.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   4.7243  4.7243  4.7243\n",
            "  ANGSTROM\n"
            "  CELL VECTORS\n" "   5.0000  0.0000  0.0000\n"
            "   0.0000  5.0000  0.0000\n" "   0.0000  0.0000  5.0000\n"
            "&END\n" "&ATOMS\n"
            "*Na.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   0.0000  0.0000  0.0000\n"
            "*Cl.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   2.5000  2.5000  2.5000\n",
            "  SCALE\n"
            "  CELL VECTORS\n" "   9.4486  0.0000  0.0000\n"
            "   0.0000  9.4486  0.0000\n" "   0.0000  0.0000  9.4486\n"
            "&END\n" "&ATOMS\n"
            "*Na.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   0.0000  0.0000  0.0000\n"
            "*Cl.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   0.5000  0.5000  0.5000\n",
            "  SCALE\n"
            "  CELL VECTORS\n" "   9.4486  0.0000  0.0000\n"
            "   0.0000  9.4486  0.0000\n" "   0.0000  0.0000  9.4486\n"
            "&END\n" "&ATOMS\n"
            "*Na.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   0.0000  0.0000  0.0000\n"
            "*Cl.uspp736.pbe BINARY\n" "  LMAX=F\n" "  1\n"
            "   0.5000  0.5000  0.5000\n"]
    tail = "ISOTOPE\n"\
           "22.99\n"\
           "35.453\n"\
           "&END\n"
    for i in range(len(body)):
        target = head + body[i] + tail
        Mol.setFmt(fmts[i])
        f = StringIO()
        cpmd.writer(Mol, f, p)
        assert f.getvalue() == target


def test_cpmd_write_kpoints_constraints():
    Mol = Molecule()
    Mol.newAtom('Na', (0, 0, 0), fix=[True])
    Mol.newAtom('Na', (1, 1, 1), fix=[True, True, True])
    p = newParam('cpmd')
    p["&ATOMS"] = "CONSTRAINTS\n"\
                  "FIX COM\n"\
                  "END CONSTRAINTS\n"
    head = "&CPMD\n"\
           "&END\n"\
           "&SYSTEM\n"\
           "  CELL VECTORS\n"\
           "   1.0000  0.0000  0.0000\n"\
           "   0.0000  1.0000  0.0000\n"\
           "   0.0000  0.0000  1.0000\n"
    tail = "&END\n"\
           "&ATOMS\n"\
           "*Na.uspp736.pbe BINARY\n"\
           "  LMAX=F\n"\
           "  2\n"\
           "   0.0000  0.0000  0.0000\n"\
           "   1.0000  1.0000  1.0000\n"\
           "CONSTRAINTS\n"\
           "FIX ATOMS\n"\
           "1\n"\
           "    2\n"\
           "FIX COORDINATES\n"\
           "1\n"\
           "1  0 1 1\n"\
           "\n"\
           "FIX COM\n"\
           "END CONSTRAINTS\n"\
           "ISOTOPE\n"\
           "22.99\n"\
           "&END\n"

    def assertKpoint():
        f = StringIO()
        cpmd.writer(Mol, f, p)
        f = f.getvalue()
        assert f == target
    Mol.setKpoints('active', 'gamma')
    target = head + tail
    assertKpoint()
    Mol.setKpoints('active', 'mpg')
    Mol.setKpoints('mpg', ('2', '0', '0', '0', '0', '0'))
    target = head +\
        "  KPOINTS MONKHORST-PACK\n  2    0    0   \n" +\
        tail
    assertKpoint()
    Mol.setKpoints('mpg', ('2', '0', '0', '0', '0', '0.5'))
    target = head +\
        "  KPOINTS MONKHORST-PACK SHIFT=0    0    0.5 \n" +\
        "  2    0    0   \n" +\
        tail
    assertKpoint()
    Mol.setKpoints('active', 'discrete')
    Mol.setKpoints('discrete',
                   [['0.0', '0.0', '0.0', '0.75'],
                    ['0.5', '0.5', '0.5', '0.25']])
    Mol.setKpoints('options', {'crystal': False, 'bands': False})
    target = head +\
        "  KPOINTS\n  2\n  0.0  0.0  0.0  0.75\n  0.5  0.5  0.5  0.25\n" +\
        tail
    assertKpoint()
    Mol.setKpoints('options', {'crystal': True, 'bands': False})
    target = head +\
        "  KPOINTS SCALED\n  2\n  0.0  0.0  0.0  0.75\n"\
        "  0.5  0.5  0.5  0.25\n" +\
        tail
    assertKpoint()
    Mol.setKpoints('discrete',
                   [['0.0', '0.0', '0.0', '2'],
                    ['0.5', '0.5', '0.5', '0']])
    Mol.setKpoints('options', {'crystal': False, 'bands': True})
    target = head +\
        "  KPOINTS BANDS\n  2    0.0  0.0  0.0  0.5  0.5  0.5 \n"\
        "  0 0. 0. 0. 0. 0. 0.\n" +\
        tail
    assertKpoint()
    Mol.setKpoints('options', {'crystal': True, 'bands': True})
    target = head +\
        "  KPOINTS SCALED BANDS\n  2    0.0  0.0  0.0  0.5  0.5  0.5 \n"\
        "  0 0. 0. 0. 0. 0. 0.\n" +\
        tail
    assertKpoint()

    Mol.newAtom('Na', (2, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (3, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (4, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (5, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (6, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (7, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (8, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (9, 1, 1), fix=[True, True, True])
    Mol.newAtom('Na', (1, 2, 1), fix=[True, True, True])
    Mol.newAtom('Na', (1, 3, 1), fix=[True, True, True])
    Mol.setKpoints('active', 'gamma')
    tail2 = "&END\n"\
            "&ATOMS\n"\
            "*Na.uspp736.pbe BINARY\n"\
            "  LMAX=F\n"\
            "  12\n"\
            "   0.0000  0.0000  0.0000\n"\
            "   1.0000  1.0000  1.0000\n"\
            "   2.0000  1.0000  1.0000\n"\
            "   3.0000  1.0000  1.0000\n"\
            "   4.0000  1.0000  1.0000\n"\
            "   5.0000  1.0000  1.0000\n"\
            "   6.0000  1.0000  1.0000\n"\
            "   7.0000  1.0000  1.0000\n"\
            "   8.0000  1.0000  1.0000\n"\
            "   9.0000  1.0000  1.0000\n"\
            "   1.0000  2.0000  1.0000\n"\
            "   1.0000  3.0000  1.0000\n"\
            "CONSTRAINTS\n"\
            "FIX ATOMS\n"\
            "11\n"\
            "    2    3    4    5    6    7    8    9   10   11\n"\
            "   12\n"\
            "FIX COORDINATES\n"\
            "1\n"\
            "1  0 1 1\n"\
            "END CONSTRAINTS\n"\
            "ISOTOPE\n"\
            "22.99\n"\
            "&END\n"
    target = head + tail2
    f = StringIO()
    cpmd.writer(Mol, f, newParam('cpmd'))
    for i in range(len(target.split('\n'))):
        print(f.getvalue().split('\n')[i], target.split('\n')[i])
        assert f.getvalue().split('\n')[i] == target.split('\n')[i]
    assert f.getvalue() == target
