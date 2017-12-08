#include <QString>
#include <QtTest>
#include <iostream>
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

bool floatComp(float a, float b)
{
    return std::abs(a-b)<std::numeric_limits<float>::epsilon();
}

void LibVipsterTest::testVec()
{
    constexpr Vec v1{1,1,1};
    constexpr Vec v2{1,2,3};
    constexpr Vec v3{1,-2,1};
    constexpr Vec v4{1.5,1.5,1.5};
    constexpr Mat m1{v1,v2,v3};
    constexpr Mat m2{Vec{1,1,1},Vec{1,2,-2},Vec{1,3,1}};
    constexpr Mat m3{Vec{4./3.,-0.5,1./6.},Vec{1./3.,0,-1./3.},Vec{-2./3.,0.5,1./6.}};
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
    QVERIFY2(v1-2 == -v1, "Vec operator- (unary)");
    QVERIFY2(floatComp(Vec_dot(v1,v2), 6), "Vec_dot");
    QVERIFY2(Vec_cross(v1,v2) == v3, "Vec_cross");
    QVERIFY2(std::abs(Vec_length(v2)-std::sqrt(14))<std::numeric_limits<float>::epsilon(), "Vec_length");
    QVERIFY2(Mat_det(m1)==6, "Mat_det");
    QVERIFY2(Mat_trans(m1)==m2, "Mat_trans");
    QVERIFY2((m1*v1==Vec{3,6,0}), "Mat operator* (vec right)");
    QVERIFY2(v1*m1==m2*v1, "Mat operator* (vec left)");
    QVERIFY2(Mat_inv(m1) == m3, "Mat_inv");
    QVERIFY_EXCEPTION_THROWN(Mat_inv({v1,v2,v4}), std::invalid_argument);
}

void LibVipsterTest::testAtom()
{
    Atom a1{"C"};
    Atom a2{"C",{{0,0,0}},0,{{false,false,false}},0};
    Atom a3{"O",{{1,2,3}}};
    QVERIFY2(a1 == a2, "Atom operator==");
    QVERIFY2(a1 != a3, "Atom operator!=");
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
    Atom atom{"C"};
    Atom atom2{"H", {0.5,0.5,0.5}, 0.5, {false, true, false}, true};
    Step step;
    /*
     * getComment,setComment
     */
    QVERIFY2(step.getComment() == "", "Step: getComment");
    step.setComment("Testäößą☭");
    QVERIFY2(step.getComment() == "Testäößą☭", "Step: getComment");
    /*
     * newAtom, getAtom, getNat
     */
    step.newAtom();
    step.newAtom({"C"});
    step.newAtom({"C", {0.,0.,0.}});
    step.newAtom({"C", {0.,0.,0.}, 0.});
    step.newAtom({"C", {0.,0.,0.}, 0., {false, false, false}});
    step.newAtom({"C", {0.,0.,0.}, 0., {false, false, false}, false});
    step.newAtom(Atom{"C",{0.,0.,0.},0.,{false,false,false},false});
    step.newAtom(atom);
    step.newAtom(atom, AtomFmt::Alat);
    QVERIFY2(step.getNat() == 9, "Step: getNat");
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "Step: newAtom (all)/ getAtom (unformatted)" + std::to_string(i);
        QVERIFY2(step.getAtom(i) == atom, msg.c_str());
    }
    /*
     * getAtoms, setAtom, delAtom, newAtoms
     */
    step.setAtom(0, atom2);
    step.setAtom(1, {"H", {0.5,0.5,0.5}, 0.5, {false, true, false}, true});
    step.setAtom(2, atom2, AtomFmt::Alat);
    for(uint i=3;i!=9;++i)
    {
        step.delAtom(3);
    }
    step.newAtoms({{atom2, atom2, atom2}});
    QVERIFY2(step.getNat() == 6 , "Step: newAtoms");
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "Step: setAtom (all)/ delAtom/ newAtoms" + std::to_string(i);
        QVERIFY2(step.getAtom(i) == atom2, msg.c_str());
    }
    for(const Atom& at:step.getAtoms())
    {
        std::string msg = "Step: getAtoms";
        QVERIFY2(at == atom2, msg.c_str());
    }
    /*
     * getCellDim, setCellDim
     */
    QVERIFY2(floatComp(step.getCellDim(), 1), "Step: getCellDim");
    step.setCellDim(2);
    QVERIFY2((step.getAtom(0).coord == Vec{0.5,0.5,0.5}), "Step: setCellDim");
    QVERIFY2(floatComp(step.getCellDim(), 2),"Step: setCellDim");
    step.setCellDim(4, true);
    QVERIFY2((step.getAtom(0).coord == Vec{1,1,1}), "Step: setCellDim (scaling)");
    QVERIFY2(floatComp(step.getCellDim(), 4),"Step: setCellDim (scaling)");
    /*
     * getCellVec, setCellVec
     */
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{1,0,0}},{{0,1,0}},{{0,0,1}} }}), "Step: getCellVec");
    step.setCellVec({{{{4,0,0}},{{0,2,0}},{{0,0,1}}}});
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{4,0,0}},{{0,2,0}},{{0,0,1}} }}), "Step: setCellVec");
    QVERIFY2((step.getAtom(0).coord == Vec{1,1,1}), "Step: setCellVec");
    step.setCellVec({{{{2,0,0}},{{0,1,0}},{{0,0,1}}}},true);
    QVERIFY2((step.getCellVec() == std::array<Vec,3>{{ {{2,0,0}},{{0,1,0}},{{0,0,1}} }}), "Step: setCellVec (scaling)");
    QVERIFY2((step.getAtom(0).coord == Vec{0.5,0.5,1}), "Step: setCellVec (scaling)");
    /*
     * formatted get/set/new
     */
    step.newAtom({"H", Vec{0.5,0.5,1}*bohrrad, 0.5, {false,true,false},true},AtomFmt::Angstrom);
    step.newAtom({"H", {0.125,0.125,0.25},0.5,{false,true,false},true}, AtomFmt::Alat);
    step.newAtom({"H", {0.0625,0.125,0.25},0.5,{false,true,false},true}, AtomFmt::Crystal);
    step.newAtoms({atom,atom,atom});
    step.setAtom(9,{"H", Vec{0.5,0.5,1}*bohrrad, 0.5, {false,true,false},true},AtomFmt::Angstrom);
    step.setAtom(10,{"H", {0.125,0.125,0.25},0.5,{false,true,false},true}, AtomFmt::Alat);
    step.setAtom(11,{"H", {0.0625,0.125,0.25},0.5,{false,true,false},true}, AtomFmt::Crystal);
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Alat).coord==Vec{0.125,0.125,0.25}), "Step: getAtom (alat)");
    for(const Atom& at:step.getAtomsFmt(AtomFmt::Alat)){
        QVERIFY2((at.coord == Vec{0.125,0.125,0.25}), "Step: getAtoms (alat)");
    }
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Crystal).coord==Vec{0.0625,0.125,0.25}), "Step: getAtom (crystal)");
    for(const Atom& at:step.getAtomsFmt(AtomFmt::Crystal)){
        QVERIFY2((at.coord == Vec{0.0625,0.125,0.25}), "Step: getAtoms (crystal)");
    }
    QVERIFY2((step.getAtomFmt(0, AtomFmt::Angstrom).coord==Vec{0.5,0.5,1}*bohrrad), "Step: getAtom (angstrom)");
    for(const Atom& at:step.getAtomsFmt(AtomFmt::Angstrom)){
        QVERIFY2((at.coord == Vec{0.5,0.5,1}*bohrrad), "Step: getAtoms (angstrom)");
    }
    step.setCellDim(4,false,AtomFmt::Angstrom);
    QVERIFY2(floatComp(step.getCellDim(), 4*invbohr), "Step: setCellDim (formatted)");
    QVERIFY2(floatComp(step.getCellDim(AtomFmt::Angstrom), 4), "Step: getCellDim (formatted)");
    step.setCellDim(4);
    /*
     * getCenter
     */
    QVERIFY2((step.getCenter(true) == Vec{0.5,0.5,1}), "Step: getCenter (of Mass)");
    QVERIFY2((step.getCenter() == Vec{4,2,2}), "Step: getCenter (of Cell)");
    /*
     * getTypes, getNtyp
     */
    QVERIFY2((step.getNtyp() == 1), "Step: getNtyp");
    step.setAtom(0,atom);
    QVERIFY2((step.getNtyp() == 2), "Step: getNtyp");
    QVERIFY2((step.getTypes() == std::set<std::string>{"H","C"}), "Step: getTypes");
}

void LibVipsterTest::testMolecule()
{
    Step step;
    Molecule mol;
    /*
     * default constructor, name, Nstep
     */
    QVERIFY2(mol.getName() == "New Molecule", "Molecule: getName");
    mol.setName("Test-…üø©");
    QVERIFY2(mol.getName() == "Test-…üø©", "Molecule: setName");
    QVERIFY2(mol.getNstep() == 1, "Molecule: getNstep");
    /*
     * non-default constructor
     */
    mol = Molecule{"Testitest", 3};
    QVERIFY2(mol.getName() == "Testitest", "Molecule: non-default constructor");
    QVERIFY2(mol.getNstep() == 3, "Molecule: non-default constructor");
    /*
     * newStep(s)
     */
    mol.newStep();
    mol.newStep(step);
    QVERIFY2(mol.getNstep() == 5, "Molecule: newStep");
    mol.newSteps({Step{},Step{},Step{}});
    QVERIFY2(mol.getNstep() == 8, "Molecule: newSteps");
    /*
     * Batch modify steps
     */
    for(Step &s: mol.getSteps()){
        QVERIFY2(floatComp(s.getCellDim(), 1), "Molecule: getSteps");
    }
    mol.setCellDimAll(2);
    for(Step& s:mol.getSteps()){
        QVERIFY2(floatComp(s.getCellDim(), 2), "Molecule: setCellDim");
    }
    /*
     * K-Points
     */
    QVERIFY2(mol.getKPoints().active == KPointFmt::Gamma, "Molecule: getKPoints");
    QVERIFY2(mol.getKPoints().mpg.x == 0, "Molecule: getKPoints");
    QVERIFY2(mol.getKPoints().discrete.properties == KPoints::Discrete::none, "Molecule: getKPoints");
    QVERIFY2(mol.getKPoints().discrete.kpoints.empty(), "Molecule: getKPoints");
    KPoints k{KPointFmt::MPG,{6,1,1,0,0,0},{KPoints::Discrete::crystal,{{Vec{1,2,3},0.5}}}};
    mol.setKPoints(k);
    QVERIFY2(mol.getKPoints().active == KPointFmt::MPG, "Molecule: setKPoints");
    mol.getKPoints().active = KPointFmt::Discrete;
    QVERIFY2(mol.getKPoints().active == KPointFmt::Discrete, "Molecule: getKPoints");
    QVERIFY2(mol.getKPoints().mpg.x == 6, "Molecule: setKPoints");
    QVERIFY2(mol.getKPoints().discrete.properties == KPoints::Discrete::crystal, "Molecule: setKPoints");
    QVERIFY2(mol.getKPoints().discrete.kpoints.size() == 1, "Molecule: setKPoints");
    QVERIFY2(floatComp(mol.getKPoints().discrete.kpoints[0].weight, 0.5), "Molecule: setKPoints");
    QVERIFY2((mol.getKPoints().discrete.kpoints[0].pos == Vec{1,2,3}), "Molecule: setKPoints");
    /*
     * const-correctness
     */
    const Molecule& mref = mol;
    mol.getStep(0).setCellDim(4);
    QVERIFY2(floatComp(mref.getStep(0).getCellDim(), 4), "Molecule: getStep const");
    mol.setCellDimAll(4);
    for(auto& s:mref.getSteps()){
        QVERIFY2(floatComp(s.getCellDim(), 4), "Molecule: getSteps const");
    }
    QVERIFY2(mol.getKPoints().active == mref.getKPoints().active, "Molecule: getKpoints const");
    QVERIFY2(mol.getKPoints().mpg.x == mref.getKPoints().mpg.x, "Molecule: getKpoints const");
    QVERIFY2(floatComp(mol.getKPoints().discrete.kpoints[0].weight, mref.getKPoints().discrete.kpoints[0].weight), "Molecule: getKpoints const");
}

QTEST_APPLESS_MAIN(LibVipsterTest)

#include "tst_libvipstertest.moc"
