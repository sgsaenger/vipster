#include <QString>
#include <QtTest>
#include <iostream>
#include <sstream>
#include <global.h>
#include <step.h>
#include <molecule.h>

using namespace Vipster;

class LibVipsterTest : public QObject
{
    Q_OBJECT

public:
    LibVipsterTest();

private Q_SLOTS:
    void testVec();
    void testAtom();
    void testPse();
    void testStep();
    void testMolecule();
};

LibVipsterTest::LibVipsterTest()
{
}

void LibVipsterTest::testVec()
{
    Vec v1{1,1,1};
    Vec v2{1,2,3};
    Vec v3{4,0,5};
    Vec v4{1.5,1.5,1.5};
    Mat m1{v1,v2,v3};
    std::ostringstream s;
    s << v4;
    QVERIFY2(v1 == v1, "Vec operator==");
    QVERIFY2(v1 != v2, "Vec operator!=");
    QVERIFY2(v1+0.5 == v4, "Vec operator+ (float right)");
    QVERIFY2(0.5+v1 == v4, "Vec operator+ (float left)");
    QVERIFY2(v4-0.5 == v1, "Vec operator- (float)");
    QVERIFY2(v1*1.5 == v4, "Vec operator* (float right)");
    QVERIFY2(1.5*v1 == v4, "Vec operator* (float left)");
    QVERIFY2(v4/1.5 == v1, "Vec operator/ (float)");
    QVERIFY2(v1+v1 == 2*v1, "Vec operator+ (Vec)");
    QVERIFY2(v4-v1 == v1/2, "Vec operator- (Vec)");
    QVERIFY2(Vec_dot(v2,v3) == 19, "Vec_dot");
    QVERIFY2(Vec_length(v2)-std::sqrt(14)<std::numeric_limits<float>::epsilon(), "Vec_length");
    QVERIFY2(s.str()=="Vec: [1.5, 1.5, 1.5]", "Vec operator<<");
    QVERIFY2(Mat_det(m1)==9, "Mat_det");
    QVERIFY_EXCEPTION_THROWN(Mat_inv({v1,v2,v4}), std::invalid_argument);
}

void LibVipsterTest::testAtom()
{

}

void LibVipsterTest::testPse()
{
    PseMap p;
    auto pseComp = [](const PseEntry& comp1, const PseEntry& comp2){
        return std::tie(comp1.bondcut, comp1.col, comp1.covr, comp1.CPNL, comp1.CPPP, comp1.m, comp1.PWPP, comp1.vdwr, comp1.Z)
                ==
               std::tie(comp2.bondcut, comp2.col, comp2.covr, comp2.CPNL, comp2.CPPP, comp2.m, comp2.PWPP, comp2.vdwr, comp2.Z);
    };

    QVERIFY2(pseComp(p["C"], Vipster::pse.at("C")), "pse direct mismatch");
    QVERIFY2(pseComp(p["c"], Vipster::pse.at("C")), "pse lower-case mismatch");
    QVERIFY2(pseComp(p["C1"], Vipster::pse.at("C")), "pse number mismatch");
    QVERIFY2(pseComp(p["CA"], Vipster::pse.at("C")), "pse letter mismatch");
    for(auto &i:p){
        QVERIFY2(pseComp(i.second, Vipster::pse.at("C")), "pse loop mismatch");
    }
    QVERIFY2(pseComp(p["Zzyzzyx"], Vipster::pse.at("X")), "pse miss mismatch");
    Molecule m;
    Step &s = m.getStep(0);
    QVERIFY2(pseComp((*m.pse)["C"], Vipster::pse.at("C")), "mol mismatch");
    QVERIFY2(pseComp((*s.pse)["C"], Vipster::pse.at("C")), "step existing mismatch");
    QVERIFY2(pseComp((*s.pse)["H"], Vipster::pse.at("H")), "step new mismatch");
}

void LibVipsterTest::testStep()
{
    Atom atom{"C", {0.,0.,0.}, 0., {false, false, false}, false};
    Atom atom2{"H", {0.5,0.5,0.5}, 0.5, {false, false, false}, false};

    Step step;

    // newAtom, getAtom, getNat
    step.newAtom();
    step.newAtom("C");
    step.newAtom("C", {0.,0.,0.});
    step.newAtom("C", {0.,0.,0.}, 0.);
    step.newAtom("C", {0.,0.,0.}, 0., {false, false, false});
    step.newAtom("C", {0.,0.,0.}, 0., {false, false, false}, false);
    step.newAtom(atom);
    step.newAtom(Atom{"C",{0.,0.,0.},0.,{false,false,false},false});
    step.newAtom(atom, AtomFmt::Alat);
    QVERIFY2(step.getNat() == 9, "step: nat mismatch");
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "step: atom mismatch at pos " + std::to_string(i);
        QVERIFY2(step.getAtom(i) == atom, msg.c_str());
    }

    // getAtoms, setAtom, delAtom
    step.setAtom(0, atom2);
    step.setAtom(1, Atom{"H",{0.5,0.5,0.5},0.5,{false,false,false},false});
    step.setAtom(2, "H", {0.5,0.5,0.5}, 0.5, {false, false, false}, false);
    step.setAtom(3, atom2, AtomFmt::Alat);
    for(uint i=4;i!=9;++i)
    {
        step.delAtom(4);
    }
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "step: atom mismatch at pos " + std::to_string(i);
        QVERIFY2(step.getAtom(i) == atom2, msg.c_str());
    }
    for(const Atom& at:step.getAtoms())
    {
        std::string msg = "step: atom mismatch";
        QVERIFY2(at == atom2, msg.c_str());
    }
    // getCellDim, setCellDim, getCellVec, setCellVec, getAtomFt, getAtomsFmt
    QVERIFY2(step.getCellDim() == 1, "step: CellDim mismatch");
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Alat).coord == Vec{0.5,0.5,0.5}), "step: Coord mismatch");
    step.setCellDim(2);
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Alat).coord == Vec{0.25,0.25,0.25}), "step: Coord mismatch");
    step.setCellDim(4, true);
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Alat).coord == Vec{0.25,0.25,0.25}), "step: Coord mismatch");
    QVERIFY2(step.getCellDim() == 4, "step: CellDim mismatch");
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{1,0,0}},{{0,1,0}},{{0,0,1}} }}), "step: CellVec mismatch");
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Crystal).coord == Vec{0.25,0.25,0.25}), "step: Coord mismatch");
    step.setCellVec(4,0,0,0,2,0,0,0,1);
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{4,0,0}},{{0,2,0}},{{0,0,1}} }}), "step: CellVec mismatch");
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Crystal).coord == Vec{0.0625,0.125,0.25}), "step: Coord mismatch");
    step.setCellVec(Vec{1,0,0},Vec{0,1,0},Vec{0,0,1},true);
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Crystal).coord == Vec{0.0625,0.125,0.25}), "step: Coord mismatch");
    for(const Atom& at:step.getAtomsFmt(AtomFmt::Crystal)){
        QVERIFY2((at.coord == Vec{0.0625,0.125,0.25}), "step: coord mismatch");
    }
    step.setCellDim(4,false,AtomFmt::Angstrom);
    QVERIFY2(step.getCellDim()!=4, "step CellDim mismatch");
    QVERIFY2(step.getCellDim(AtomFmt::Angstrom)==4, "step CellDim mismatch");
    step.setCellDim(4);
    QVERIFY2(step.getCellDim()==4, "step CellDim mismatch");
    // getCenter
//    QVERIFY2((step.getCenter(true) == Vec{0.125,0.25,0.5}), "step: Center mismatch");
    QVERIFY2((step.getCenter() == Vec{2,2,2}), "step: Center mismatch");
    // getTypes, getNtyp
    QVERIFY2((step.getNtyp() == 1), "step: Ntyp mismatch");
    step.setAtom(0);
    QVERIFY2((step.getNtyp() == 2), "step: Ntyp mismatch");
    QVERIFY2((step.getTypes() == std::set<std::string>{"H","C"}), "step: types mismatch");
}

void LibVipsterTest::testMolecule()
{
    Molecule mol1;
    QVERIFY2(mol1.getName() == "New Molecule", "mol1: name mismatch");
    QVERIFY2(mol1.getNstep() == 1, "mol1: length mismatch");
    Molecule mol2{"Test molecule"};
    QVERIFY2(mol2.getName() == "Test molecule", "mol2: name mismatch");
    QVERIFY2(mol2.getNstep() == 1, "mol2: length mismatch");
    Molecule mol3{"Test mol3", 3};
    QVERIFY2(mol3.getName() == "Test mol3", "mol3: name mismatch");
    QVERIFY2(mol3.getNstep() == 3, "mol3: length mismatch");
    for(Step &s: mol3.getSteps()){
        QVERIFY2(s.getCellDim()==1, "mol3: CellDim mismatch");
    }
    mol3.setCellDimAll(2);
    for(Step& s:mol3.getSteps()){
        QVERIFY2(s.getCellDim()==2, "mol3: CellDim mismatch");
    }
    Molecule mol4 = Molecule("Move test");
    QVERIFY2(mol4.getName() == "Move test", "mol4: name mismatch");
    QVERIFY2(mol4.getNstep() == 1, "mol4: length mismatch");
}

QTEST_APPLESS_MAIN(LibVipsterTest)

#include "tst_libvipstertest.moc"
