from test_preamble import *
from vipster import Molecule


def test_atoms_bohr():
    Mol = Molecule()
    Mol.setCellDim(2)
    Mol.newAtom()
    Mol.newAtom('H', (1, 1, 1), 1, [True, False], True)
    assert Mol.nat == 2
    assert Mol.ntyp == 2
    assert Mol.getTypes() == ['H', 'C']
    assert Mol.getCellDim() == 2
    assert atom_equal(Mol.getAtom(0, charge=True, fix=True, hidden=True),
                      ['C', (0, 0, 0), 0., [False, False, False], False])
    assert atom_equal(Mol.getAtom(1, charge=True, fix=True, hidden=True),
                      ['H', (1, 1, 1), 1., [True, False, False], True])
    assert list(map(len, Mol.getBonds(1.1))) == [1, 1, 1, 1, 1, 1, 1, 1]
    Mol.delAtom(1)
    assert Mol.nat == 1
    assert Mol.getTypes() == ['C']
    Mol.newAtoms(2)
    assert Mol.getTypes() == ['X', 'C']
    Mol.setAtom(1, 'O', (0.5, 0.5, 0.5), 2, [False, True], False)
    Mol.setAtom(2, 'U', (1, 1, 1), 2, [True], False)
    assert atom_equal(Mol.getAtom(1, charge=True, fix=True, hidden=True),
                      ['O', (0.5, 0.5, 0.5), 2., [False, True, False], False])
    assert atom_equal(Mol.getAtom(2, charge=True, fix=True, hidden=True),
                      ['U', (1, 1, 1), 2., [True, False, False], False])
    assert Mol.getTypes() == ['C', 'O', 'U']
    target = (['C', (0, 0, 0), 0., [False, False, False], False],
              ['O', (0.5, 0.5, 0.5), 2., [False, True, False], False],
              ['U', (1, 1, 1), 2., [True, False, False], False])
    for a, t in zip(Mol.getAtoms(True, True, True), target):
        assert atom_equal(a, t)


def test_atoms_fmt():
    Mol = Molecule()
    # standard: bohr:
    assert Mol.getFmt() == 'bohr'
    assert float_equal(Mol.getCellDim(), 1)
    assert vec_equal(Mol.getVec(), ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    Mol.newAtom('C', (1, 1, 1))
    assert atom_equal(Mol.getAtom(0), ('C', (1, 1, 1)))
    # transition to angstrom:
    Mol.setFmt("angstrom", scale=True)
    assert Mol.getFmt() == 'angstrom'
    Mol.setCellDim(1)
    assert float_equal(Mol.getCellDim(), 1)
    assert atom_equal(Mol.getAtom(0, fmt="angstrom"), ('C', (1, 1, 1)))
    Mol.setCellDim(2.5, scale=True, fmt='angstrom')
    Mol.newAtom('C', (2.5, 2.5, 2.5), fmt='angstrom')
    for at in Mol.getAtoms():
        assert atom_equal(at, ('C', (2.5, 2.5, 2.5)))
    Mol.setVec(((2, 0, 0), (0, 2, 0), (0, 0, 2)), scale=True)
    assert vec_equal(Mol.getVec(), ((2, 0, 0), (0, 2, 0), (0, 0, 2)))
    Mol.newAtom('C', (5, 5, 5), fmt='angstrom')
    Mol.newAtom('C', (2, 2, 2), fmt='alat')
    Mol.newAtom('C', (1, 1, 1), fmt='crystal')
    for at in Mol.getAtoms():
        assert atom_equal(at, ('C', (5, 5, 5)))
    for at in Mol.getAtoms(fmt='alat'):
        assert atom_equal(at, ('C', (2, 2, 2)))
    for at in Mol.getAtoms(fmt='crystal'):
        assert atom_equal(at, ('C', (1, 1, 1)))
    # scale back down to bohr
    Mol.setCellDim(2.5, scale=True, fmt='bohr')
    for at in Mol.getAtoms(fmt='bohr'):
        assert atom_equal(at, ('C', (5, 5, 5)))
    Mol.setFmt('alat', scale=True)
    for at in Mol.getAtoms(fmt='bohr'):
        assert atom_equal(at, ('C', (12.5, 12.5, 12.5)))
    Mol.setFmt('crystal', scale=True)
    for at in Mol.getAtoms(fmt='bohr'):
        assert atom_equal(at, ('C', (62.5, 62.5, 62.5)))
    assert vec_equal(Mol.getCenter(), (2.5, 2.5, 2.5))
    assert vec_equal(Mol.getCenter(com=True), (62.5, 62.5, 62.5))


def test_modify():
    Mol = Molecule()
    Mol.setVec(((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    Mol.newAtom('C', (0.5, 0.5, 0.5), fmt='crystal')
    assert vec_equal(Mol.getVec(), ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    assert Mol.nat == 1
    Mol.mult(2, 1, 1)
    assert vec_equal(Mol.getVec(), ((10, 0, 0), (0, 5, 0), (0, 0, 5)))
    assert Mol.nat == 2
    Mol.mult(1, 2, 2)
    assert vec_equal(Mol.getVec(), ((10, 0, 0), (0, 10, 0), (0, 0, 10)))
    assert Mol.nat == 8
    Mol.setVec(((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    Mol.crop()
    assert Mol.nat == 1
    Mol.mult(2, 2, 2)
    assert vec_equal(Mol.getVec(), ((10, 0, 0), (0, 10, 0), (0, 0, 10)))
    assert Mol.nat == 8
    Mol.setVec(((7, 0, 0), (0, 7, 0), (0, 0, 7)))
    Mol.wrap()
    assert Mol.nat == 8
    assert list(map(len, Mol.getBonds(1.1))) == [24, 0, 0, 0, 0, 0, 0, 0]
    assert Mol.getUndo() == 'wrap atoms'
    Mol.undo()
    assert list(map(len, Mol.getBonds(1.1))) == [0, 4, 4, 4, 4, 4, 4, 0]
    Mol.wrap()
    Mol.reshape(((14, 7, 0), (0, 7, 0), (0, 0, 7)))
    assert list(map(len, Mol.getBonds(1.1))) == [38, 0, 10, 0, 0, 0, 0, 0]
    assert vec_equal(Mol.getVec(), ((14, 7, 0), (0, 7, 0), (0, 0, 7)))
    Mol.undo()
    Mol.setVec(((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    assert atom_equal(Mol.getAtom(1, fmt="crystal"), ('C', (0.1, 0.5, 0.5)))
    Mol.align(1, 'x')
    assert vec_equal(Mol.getVec(), ((0, -5, 0), (5, 0, 0), (0, 0, 5)))
    Mol.undo()
    Mol.align(1, 'y')
    assert vec_equal(Mol.getVec(), ((5, 0, 0), (0, 5, 0), (0, 0, 5)))
    Mol.align(1, 'z')
    assert vec_equal(Mol.getVec(), ((5, 0, 0), (0, 0, 5), (0, -5, 0)))
    Mol.undo()
    try:
        Mol.align(0, 1)
    except ValueError:
        assert True
    except:
        assert False, "Wrong error"
    else:
        assert False, "Error missing"


def test_selection():
    Mol = Molecule()
    Mol.newAtom('C', (0, 0, 0))
    Mol.addSelection([1, (0, 0, 0)])
    assert Mol.getSelection() == [[1, (0, 0, 0)]]
    Mol.delSelection()
    assert Mol.getSelection() == []


def test_scripting():
    Mol = Molecule()
    Mol.setFmt('angstrom')
    Mol.newAtom('C', (0, 0, 0))
    Mol.newAtom('O', (1, 1, 0))
    Mol.newAtom('B', (2, 2, 0))
    Mol.newAtom('U', (0, 2, 0))
    Mol.evalScript("def g all tU\n"
                   "shi g 1-3 2")
    assert atom_equal(Mol.getAtom(3), ('U', (2, 0, 0)))
    Mol.undo()
    Mol.evalScript("rot 3 90 2-0 1")
    assert atom_equal(Mol.getAtom(3), ('U', (1, 1, np.sqrt(2))))
    Mol.undo()
    Mol.evalScript("rot 3 90 2-0\n"
                   "mir 3 (1,0,0) (0,0,1,'angstrom') (1,0,0)")
    assert atom_equal(Mol.getAtom(3), ('U', (1, -1, np.sqrt(2))))
    Mol = Molecule()
    Mol.setFmt('angstrom')
    Mol.setCellDim(15)
    Mol.newAtom('C1', (0, 0, 0))
    Mol.newAtom('C1', (1, 0, 0))
    Mol.newAtom('C1', (1, 0, 0))
    Mol.newAtom('C2', (0, 0, 3))
    Mol.newAtom('C2', (1, 0, 3))
    Mol.newAtom('C2', (1, 0, 3))
    Mol.evalScript("rot [2, 5] 60 (0, 0, 1)")
    Mol.evalScript("rot 3-5 90 (1, 0, 0) 3")
    assert len(Mol.getBonds(1.1)[0]) == 6
    Mol.evalScript("sel all tC2\npar sel sel 0-2")
    atoms = Mol.getAtoms()
    assert vec_equal(atoms[0][1], atoms[3][1]-(0,0,3))
    assert vec_equal(atoms[1][1], atoms[4][1]-(0,0,3))
    assert vec_equal(atoms[2][1], atoms[5][1]-(0,0,3))
    Mol.evalScript("rot sel 45 (1, 0, 0) 3")
    Mol.evalScript("pshift 0-2 (0,0,1) sel")
    assert atom_equal(Mol.getAtom(0), ['C1', (0, np.sqrt(2)/2, -np.sqrt(2)/2)])
    assert atom_equal(Mol.getAtom(1), ['C1', (1, np.sqrt(2)/2, -np.sqrt(2)/2)])
    assert atom_equal(Mol.getAtom(2), ['C1', (0.5, (np.sqrt(3)+np.sqrt(2))/2, -np.sqrt(2)/2)])
