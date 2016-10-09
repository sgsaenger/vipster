#include <QString>
#include <QtTest>
#include "libvipster.h"
#include <iostream>

using namespace Vipster;

class LibVipsterTest : public QObject
{
    Q_OBJECT

public:
    LibVipsterTest();

private Q_SLOTS:
    void testMolecule();
    void testStep();
};

LibVipsterTest::LibVipsterTest()
{
}

void LibVipsterTest::testMolecule()
{
    Molecule mol1{"Test molecule"};
    Molecule mol3{"Test mol3", 3};
    QVERIFY2(mol1.name == "Test molecule", "mol1: name mismatch");
    QVERIFY2(mol1.stepIdx == 0, "mol1: stepIdx mismatch");
    QVERIFY2(mol1.steps.size() == 1, "mol1: length mismatch");
    QVERIFY2(mol3.name == "Test mol3", "mol3: name mismatch");
    QVERIFY2(mol3.stepIdx == 2, "mol3: stepIdx mismatch");
    QVERIFY2(mol3.steps.size() == 3, "mol3: length mismatch");
    QVERIFY2(&mol3.curStep() == &mol3.steps[2], "mol3: step mismatch");
    for(Step& s:mol3.steps){
        QVERIFY2(s.getCellDim()==1, "mol3: CellDim mismatch");
    }
    mol3.setCellDimAll(2);
    for(Step& s:mol3.steps){
        QVERIFY2(s.getCellDim()==2, "mol3: CellDim mismatch");
    }
}

void LibVipsterTest::testStep()
{
    Molecule mol{"Test Molecule", 2};
    QVERIFY2(mol.steps.size() == 2, "mol: length mismatch");
    QVERIFY2(mol.stepIdx == 1, "mol: stepIdx mismatch");
    Step step = mol.curStep();
    Atom atom{"C", {0.,0.,0.}, 0., {false, false, false}, false};
    Atom atom2{"H", {0.5,0.5,0.5}, 0.5, {false, false, false}, false};
    auto atomComp = [atom](const Atom& comp1, const Atom& comp2){
        return std::tie(comp1.name, comp1.coord, comp1.charge, comp1.fix, comp1.hidden)
                ==
               std::tie(comp2.name, comp2.coord, comp2.charge, comp2.fix, comp2.hidden);
    };
    // newAtom, getAtom, getNat
    step.newAtom();
    step.newAtom("C");
    step.newAtom("C", {0.,0.,0.});
    step.newAtom("C", {0.,0.,0.}, 0.);
    step.newAtom("C", {0.,0.,0.}, 0., {false, false, false});
    step.newAtom("C", {0.,0.,0.}, 0., {false, false, false}, false);
    step.newAtom(atom);
    step.newAtom(Atom{"C",{0.,0.,0.},0.,{false,false,false},false});
    step.newAtom(atom, Fmt::Alat);
    QVERIFY2(step.getNat() == 9, "step: nat mismatch");
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "step: atom mismatch at pos " + std::to_string(i);
        QVERIFY2(atomComp(step.getAtom(i), atom), msg.c_str());
    }
    mol.stepIdx = 0;
    // newAtoms, getAtoms, setAtom, delAtom
    step = mol.curStep();
    QVERIFY2(step.getNat() == 0, "step: nat mismatch");
    step.newAtoms(5);
    QVERIFY2(step.getNat() == 5, "step: nat mismatch");
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "step: atom mismatch at pos " + std::to_string(i);
        QVERIFY2(atomComp(step.getAtom(i), atom), msg.c_str());
    }
    step.setAtom(0, atom2);
    step.setAtom(1, Atom{"H",{0.5,0.5,0.5},0.5,{false,false,false},false});
    step.setAtom(2, "H", {0.5,0.5,0.5}, 0.5, {false, false, false}, false);
    step.setAtom(3, atom2, Fmt::Alat);
    step.delAtom(4);
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "step: atom mismatch at pos " + std::to_string(i);
        QVERIFY2(atomComp(step.getAtom(i), atom2), msg.c_str());
    }
    for(const Atom& at:step.getAtoms())
    {
        std::string msg = "step: atom mismatch";
        QVERIFY2(atomComp(at, atom2), msg.c_str());
    }
    // getCellDim, setCellDim, getCellVec, setCellVec, getAtomFt, getAtomsFmt
    QVERIFY2(step.getCellDim() == 1, "step: CellDim mismatch");
    QVERIFY2((step.getAtomFmt(0, Fmt::Alat).coord == Vec{0.5,0.5,0.5}), "step: Coord mismatch");
    step.setCellDim(2);
    QVERIFY2((step.getAtomFmt(0, Fmt::Alat).coord == Vec{0.25,0.25,0.25}), "step: Coord mismatch");
    step.setCellDim(4, true);
    QVERIFY2((step.getAtomFmt(0, Fmt::Alat).coord == Vec{0.25,0.25,0.25}), "step: Coord mismatch");
    QVERIFY2(step.getCellDim() == 4, "step: CellDim mismatch");
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{1,0,0}},{{0,1,0}},{{0,0,1}} }}), "step: CellVec mismatch");
    QVERIFY2((step.getAtomFmt(0, Fmt::Crystal).coord == Vec{0.25,0.25,0.25}), "step: Coord mismatch");
    step.setCellVec(4,0,0,0,2,0,0,0,1);
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{4,0,0}},{{0,2,0}},{{0,0,1}} }}), "step: CellVec mismatch");
    QVERIFY2((step.getAtomFmt(0, Fmt::Crystal).coord == Vec{0.0625,0.125,0.25}), "step: Coord mismatch");
    step.setCellVec(Vec{1,0,0},Vec{0,1,0},Vec{0,0,1},true);
    QVERIFY2((step.getAtomFmt(0, Fmt::Crystal).coord == Vec{0.0625,0.125,0.25}), "step: Coord mismatch");
    for(const Atom& at:step.getAtomsFmt(Fmt::Crystal)){
        QVERIFY2((at.coord == Vec{0.0625,0.125,0.25}), "step: coord mismatch");
    }
    step.setCellDim(4,false,Fmt::Angstrom);
    QVERIFY2(step.getCellDim()!=4, "step CellDim mismatch");
    QVERIFY2(step.getCellDim(Fmt::Angstrom)==4, "step CellDim mismatch");
    step.setCellDim(4);
    QVERIFY2(step.getCellDim()==4, "step CellDim mismatch");
    // getCenter
    QVERIFY2((step.getCenter(true) == Vec{0.125,0.25,0.5}), "step: Center mismatch");
    QVERIFY2((step.getCenter() == Vec{2,2,2}), "step: Center mismatch");
    // getTypes, getNtyp
    QVERIFY2((step.getNtyp() == 1), "step: Ntyp mismatch");
    step.setAtom(0);
    QVERIFY2((step.getNtyp() == 2), "step: Ntyp mismatch");
    QVERIFY2((step.getTypes() == std::set<std::string>{"H","C"}), "step: types mismatch");
}

QTEST_APPLESS_MAIN(LibVipsterTest)

#include "tst_libvipstertest.moc"
