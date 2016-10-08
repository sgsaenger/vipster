#include <QString>
#include <QtTest>
#include "libvipster.h"

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
    Vipster::Molecule mol1{"Test molecule"};
    Vipster::Molecule mol3{"Test mol3", 3};
    QVERIFY2(mol1.name == "Test molecule", "mol1: name mismatch");
    QVERIFY2(mol1.stepIdx == 0, "mol1: stepIdx mismatch");
    QVERIFY2(mol1.steps.size() == 1, "mol1: length mismatch");
    QVERIFY2(mol3.name == "Test mol3", "mol3: name mismatch");
    QVERIFY2(mol3.stepIdx == 2, "mol3: stepIdx mismatch");
    QVERIFY2(mol3.steps.size() == 3, "mol3: length mismatch");
    QVERIFY2(&mol3.curStep() == &mol3.steps[2], "mol3: step mismatch");
}

void LibVipsterTest::testStep()
{
    Vipster::Molecule mol{"Test Molecule", 2};
    QVERIFY2(mol.steps.size() == 2, "mol: length mismatch");
    Vipster::Step step = mol.curStep();
    Vipster::Atom atom{"C", {0.,0.,0.}, 0., {false, false, false}, false};
    auto atomComp = [atom](Vipster::Atom comp){
        return std::tie(atom.name, atom.coord, atom.charge, atom.fix, atom.hidden)
                ==
               std::tie(comp.name, comp.coord, comp.charge, comp.fix, comp.hidden);
    };
    step.newAtom();
    step.newAtom("C");
    step.newAtom("C", {0.,0.,0.});
    step.newAtom("C", {0.,0.,0.}, 0.);
    step.newAtom("C", {0.,0.,0.}, 0., {false, false, false});
    step.newAtom("C", {0.,0.,0.}, 0., {false, false, false}, false);
    step.newAtom(atom);
    step.newAtom(Vipster::Atom{"C",{0.,0.,0.},0.,{false,false,false},false});
    QVERIFY2(step.getNat() == 8, "step: nat mismatch");
    for(uint i=0;i!=step.getNat();++i)
    {
        std::string msg = "step: atom mismatch at pos " + std::to_string(i);
        QVERIFY2(atomComp(step.getAtom(i)), msg.c_str());
    }
}

QTEST_APPLESS_MAIN(LibVipsterTest)

#include "tst_libvipstertest.moc"
