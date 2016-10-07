from vipster.iowrapper import _paramdict
from vipster import *
from test_preamble import *
from os import remove


def test_params():
    assert availParam() == _paramdict.keys()
    assert availParam("Fail") == None
    for i in _paramdict.keys():
        assert availParam(i) == _paramdict[i].keys()
        assert newParam(i) == _paramdict[i]['default']


def test_basic_io():
    Mol = Molecule()
    Mol.setFmt('angstrom')
    Mol.newAtom('C', (0, 0, 0))
    Mol.newAtom('O', (2.5, 2.5, 2.5))
    writeFile(Mol, 'xyz', 'test.xyz')
    Mol2,_ = readFile('test.xyz', 'xyz')
    remove('test.xyz')
    assert Mol.nat == Mol2.nat
    assert Mol.ntyp == Mol2.ntyp
    assert Mol.getTypes() == Mol2.getTypes()
    for i in range(Mol.nat):
        assert atom_equal(Mol.getAtom(i), Mol2.getAtom(i))


def test_io_errorhandling():
    Mol = Molecule()
    Mol.setFmt('angstrom')
    Mol.newAtom('C', (0, 0, 0))
    Mol.newAtom('O', (2.5, 2.5, 2.5))
    p = newParam('pwi')
    p['&control']['calculation'] = "'relax'"

    def assertFail(fmt, path, ex, parm=None):
        try:
            writeFile(Mol, fmt, path, parm)
        except ex:
            assert True
        except:
            assert False, "unexpected exception"
        else:
            assert False, "missing exception"
    assertFail('xyz', '/text.xyz', PermissionError)
    assertFail('pwi', 'test.pwi', TypeError)
    assertFail('pwi', 'test.pwi', KeyError, p)
