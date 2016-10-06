from vipster.ioplugins import lammpsData
from vipster import Molecule, newParam
from io import StringIO
from test_preamble import *


def test_lmp_parse_atomstyles():
    name = "lmp_parse"
    head = ["\n",
            "4 atoms\n",
            "2 atom types\n",
            "\n",
            "0 5 xlo xhi\n",
            "0 5 ylo yhi\n",
            "0 5 zlo zhi\n",
            "\n",
            "Masses\n",
            "\n",
            "1 22.99 #Na\n",
            "2 34.99 #Cl\n",
            "\n"]
    atoms = [["Atoms # atomic\n",
              "\n",
              "1 1 0.0 0.0 0.0\n",
              "2 1 2.5 2.5 0.0\n",
              "3 2 0.0 2.5 0.0\n",
              "4 2 2.5 0.0 0.0\n",
              "\n"],
             ["Atoms\n",
              "\n",
              "1 1 0.0 0.0 0.0\n",
              "2 1 2.5 2.5 0.0\n",
              "3 2 0.0 2.5 0.0\n",
              "4 2 2.5 0.0 0.0\n",
              "\n"],
             ["Atoms\n",
              "\n",
              "1 1 1 0.0 0.0 0.0\n",
              "2 1 1 2.5 2.5 0.0\n",
              "3 1 2 0.0 2.5 0.0\n",
              "4 1 2 2.5 0.0 0.0\n",
              "\n"],
             ["Atoms # bond\n",
              "\n",
              "1 1 1 0.0 0.0 0.0\n",
              "2 1 1 2.5 2.5 0.0\n",
              "3 1 2 0.0 2.5 0.0\n",
              "4 1 2 2.5 0.0 0.0\n",
              "\n"],
             ["Atoms\n",
              "\n",
              "1 1 1  1 0.0 0.0 0.0\n",
              "2 1 1  1 2.5 2.5 0.0\n",
              "3 1 2 -1 0.0 2.5 0.0\n",
              "4 1 2 -1 2.5 0.0 0.0\n",
              "\n"],
             ["Atoms # full\n",
              "\n",
              "1 1 1  1 0.0 0.0 0.0\n",
              "2 1 1  1 2.5 2.5 0.0\n",
              "3 1 2 -1 0.0 2.5 0.0\n",
              "4 1 2 -1 2.5 0.0 0.0\n",
              "\n"],
             ]
    for a in atoms:
        Mol, _ = lammpsData.parser(name, head + a)
        assert Mol.nat == 4
        assert Mol.ntyp == 2
        assert Mol.getTypes() == ['Na', 'Cl']
        assert Mol.pse['Cl']['m'] == 34.99
        assert float_equal(Mol.getCellDim(fmt='angstrom'), 1)
        assert vec_equal(Mol.getVec(), ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
        assert list(map(len, Mol.getBonds(1.1))) == [4, 2, 2, 0, 0, 0, 0, 0]
        assert atom_equal(Mol.getAtom(0), ['Na', [0, 0, 0]])
        assert atom_equal(Mol.getAtom(1), ['Na', [2.5, 2.5, 0]])
        assert atom_equal(Mol.getAtom(2), ['Cl', [0, 2.5, 0]])
        assert atom_equal(Mol.getAtom(3), ['Cl', [2.5, 0, 0]])
    assert atom_equal(Mol.getAtom(0, charge=True), ['Na', [0, 0, 0], 1.])
    assert atom_equal(Mol.getAtom(2, charge=True), ['Cl', [0, 2.5, 0], -1.0])


def test_lmp_parse_image_skew_fail():
    name = 'lmp_parse_other'
    target = ["\n",
              "4 atoms\n",
              "2 atom types\n",
              "\n",
              "0 5 xlo xhi\n",
              "0 5 ylo yhi\n",
              "0 5 zlo zhi\n",
              "1 2 3 xy xz yz\n",
              "\n",
              "Masses\n",
              "\n",
              "1 22.99 #Na\n",
              "2 34.99 #Cl\n",
              "\n",
              "Atoms\n",
              "\n",
              "1 1 0.0 0.0 0.0 0 0 0\n",
              "2 1 2.5 2.5 0.0 0 0 0\n",
              "3 2 0.0 2.5 0.0 0 0 0\n",
              "4 2 2.5 0.0 0.0 0 0 0\n",
              "\n"]
    Mol, _ = lammpsData.parser(name, target)
    assert Mol.nat == 4
    assert float_equal(Mol.getCellDim(fmt='angstrom'), 1)
    assert vec_equal(Mol.getVec(), ((5, 0, 0), (1, 5, 0), (2, 3, 5)))
    failtarget = ["\n",
                  "4 atoms\n",
                  "2 atom types\n",
                  "\n",
                  "0 5 xlo xhi\n",
                  "0 5 ylo yhi\n",
                  "0 5 zlo zhi\n",
                  "1 2 3 xy xz yz\n",
                  "\n",
                  "Masses\n",
                  "\n",
                  "1 22.99\n",
                  "2 34.99\n",
                  "\n",
                  "Atoms\n",
                  "\n",
                  "1 1 0.0 0.0 0.0 0 0 0\n",
                  "2 1 2.5 2.5 0.0 0 0 0\n",
                  "3 2 0.0 2.5 0.0 0 0 0\n",
                  "4 2 2.5 0.0 0.0 0 0 0\n",
                  "\n"]
    try:
        Mol, _ = lammpsData.parser(name, failtarget)
    except NotImplementedError:
        assert True
    except:
        assert False
    else:
        assert False


def test_lmp_write_proper():
    Mol = Molecule()
    Mol.setFmt('angstrom')
    Mol.setCellDim(1)
    Mol.setVec(((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    Mol.newAtom('Na', (0.0, 0.0, 0.0), charge=1.)
    Mol.newAtom('Na', (2.5, 2.5, 0.0), charge=1.)
    Mol.newAtom('Na', (0.0, 2.5, 2.5), charge=1.)
    Mol.newAtom('Na', (2.5, 0.0, 2.5), charge=1.)
    Mol.newAtom('Cl', (0.0, 2.5, 0.0), charge=-1.)
    Mol.newAtom('Cl', (2.5, 0.0, 0.0), charge=-1.)
    Mol.newAtom('Cl', (0.0, 0.0, 2.5), charge=-1.)
    Mol.newAtom('Cl', (2.5, 2.5, 2.5), charge=-1.)
    p = newParam('lmp')
    p['atom_style'] = 'full'
    p['bonds'] = True
    p['angles'] = True
    p['dihedrals'] = True
    p['impropers'] = True
    head = "\n"\
           "8 atoms\n"\
           "2 atom types\n"\
           "12 bonds\n"\
           "1 bond types\n"\
           "#1 Cl Na\n"\
           "24 angles\n"\
           "2 angle types\n"\
           "#1 Cl Na Cl\n"\
           "#2 Na Cl Na\n"\
           "48 dihedrals\n"\
           "1 dihedral types\n"\
           "#1 Cl Na Cl Na\n"\
           "8 impropers\n"\
           "2 improper types\n"\
           "#1 Cl Na Na Na\n"\
           "#2 Na Cl Cl Cl\n"\
           "\n"\
           " 0.0000  5.0000 xlo xhi\n"\
           " 0.0000  5.0000 ylo yhi\n"\
           " 0.0000  5.0000 zlo zhi\n"\
           "\n"\
           "Masses\n"\
           "\n"\
           "1 22.9900 #Na\n"\
           "2 35.4530 #Cl\n"\
           "\n"
    tail = "\nBonds\n\n"\
           "1 1 1 5\n2 1 1 6\n3 1 1 7\n4 1 2 5\n5 1 2 6\n"\
           "6 1 2 8\n7 1 3 5\n8 1 3 7\n9 1 3 8\n10 1 4 6\n"\
           "11 1 4 7\n12 1 4 8\n"\
           "\nAngles\n\n"\
           "1 1 5 1 6\n2 1 5 1 7\n3 2 1 5 2\n4 2 1 5 3\n"\
           "5 1 6 1 7\n6 2 1 6 2\n7 2 1 6 4\n8 2 1 7 3\n"\
           "9 2 1 7 4\n10 1 5 2 6\n11 1 5 2 8\n12 2 2 5 3\n"\
           "13 1 6 2 8\n14 2 2 6 4\n15 2 2 8 3\n16 2 2 8 4\n"\
           "17 1 5 3 7\n18 1 5 3 8\n19 1 7 3 8\n20 2 3 7 4\n"\
           "21 2 3 8 4\n22 1 6 4 7\n23 1 6 4 8\n24 1 7 4 8\n"\
           "\nDihedrals\n\n"\
           "1 1 6 1 5 2\n2 1 6 1 5 3\n3 1 5 1 6 2\n4 1 5 1 6 4\n"\
           "5 1 7 1 5 2\n6 1 7 1 5 3\n7 1 5 1 7 3\n8 1 5 1 7 4\n"\
           "9 1 1 5 2 6\n10 1 1 5 2 8\n11 1 1 5 3 7\n12 1 1 5 3 8\n"\
           "13 1 7 1 6 2\n14 1 7 1 6 4\n15 1 6 1 7 3\n16 1 6 1 7 4\n"\
           "17 1 1 6 2 5\n18 1 1 6 2 8\n19 1 1 6 4 7\n20 1 1 6 4 8\n"\
           "21 1 1 7 3 5\n22 1 1 7 3 8\n23 1 1 7 4 6\n24 1 1 7 4 8\n"\
           "25 1 6 2 5 3\n26 1 5 2 6 4\n27 1 8 2 5 3\n28 1 5 2 8 3\n"\
           "29 1 5 2 8 4\n30 1 2 5 3 7\n31 1 2 5 3 8\n32 1 8 2 6 4\n"\
           "33 1 6 2 8 3\n34 1 6 2 8 4\n35 1 2 6 4 7\n36 1 2 6 4 8\n"\
           "37 1 2 8 3 5\n38 1 2 8 3 7\n39 1 2 8 4 6\n40 1 2 8 4 7\n"\
           "41 1 5 3 7 4\n42 1 5 3 8 4\n43 1 8 3 7 4\n44 1 7 3 8 4\n"\
           "45 1 3 7 4 6\n46 1 3 7 4 8\n47 1 3 8 4 6\n48 1 3 8 4 7\n"\
           "\nImpropers\n\n"\
           "1 2 1 5 6 7\n2 2 2 5 6 8\n3 2 3 5 7 8\n4 2 4 6 7 8\n"\
           "5 1 5 1 2 3\n6 1 6 1 2 4\n7 1 7 1 3 4\n8 1 8 2 3 4\n"\
           "\n"
    atoms = "Atoms # full\n"\
            "\n"\
            "1 1 1  1.0000  0.0000  0.0000  0.0000\n"\
            "2 1 1  1.0000  2.5000  2.5000  0.0000\n"\
            "3 1 1  1.0000  0.0000  2.5000  2.5000\n"\
            "4 1 1  1.0000  2.5000  0.0000  2.5000\n"\
            "5 1 2 -1.0000  0.0000  2.5000  0.0000\n"\
            "6 1 2 -1.0000  2.5000  0.0000  0.0000\n"\
            "7 1 2 -1.0000  0.0000  0.0000  2.5000\n"\
            "8 1 2 -1.0000  2.5000  2.5000  2.5000\n"
    f = StringIO()
    lammpsData.writer(Mol, f, p)
    target = head + atoms + tail
    assert f.getvalue() == target


def test_lmp_write_skew():
    Mol = Molecule()
    Mol.setFmt('angstrom')
    Mol.setCellDim(1)
    Mol.setVec(((5, 0, 0), (1, 5, 0), (2, 3, 5)))
    Mol.newAtom('Na', (0, 0, 0))
    p = newParam('lmp')
    f = StringIO()
    lammpsData.writer(Mol, f, p)
    print(f.getvalue())
    target = "\n"\
             "1 atoms\n"\
             "1 atom types\n"\
             "\n"\
             " 0.0000  5.0000 xlo xhi\n"\
             " 0.0000  5.0000 ylo yhi\n"\
             " 0.0000  5.0000 zlo zhi\n"\
             " 1.0000  2.0000  3.0000 xy xz yz\n"\
             "\n"\
             "Masses\n"\
             "\n"\
             "1 22.9900 #Na\n"\
             "\n"\
             "Atoms # atomic\n"\
             "\n"\
             "1 1  0.0000  0.0000  0.0000\n"\
             "\n"
    assert f.getvalue() == target
    Mol.setVec(((5, 3, 2), (0, 5, 1), (0, 0, 5)))
    try:
        lammpsData.writer(Mol, f, p)
    except ValueError:
        assert True
    except:
        assert False
    else:
        assert False
