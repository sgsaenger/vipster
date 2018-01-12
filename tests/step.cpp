#include "atomproper.h"
#include "stepproper.h"
#include <iostream>

using namespace Vipster;

static std::ostream& operator<<(std::ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static std::ostream& operator<<(std::ostream& out, const FixVec& v)
{
    out << "FixVec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static std::ostream& operator<<(std::ostream& os, const Atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.charge
       << ", " << at.fix << ", " << at.hidden << "}";
    return os;
}

#include "catch.hpp"

TEST_CASE( "Vipster::Atom", "[atom]" ) {
    AtomProper ap{};

    SECTION( "Atom::Propref" ) {
        SECTION("Atom::name: wrapped std::string") {
            bool mod{};
            std::string s{};
            Atom::PropRef<std::string> pr = {&s, &mod};
            REQUIRE(pr == "");
            REQUIRE(mod == false);
            pr = "test";
            REQUIRE(s == "test");
            REQUIRE(s.c_str() == pr.c_str());
            REQUIRE(mod == true);
        }
        SECTION("Atom::coord: wrapped Vipster::Vec") {
            bool mod{};
            Vec v{};
            Atom::PropRef<Vec> pr = {&v, &mod};
            REQUIRE(pr == Vec{});
            REQUIRE(mod == false);
            pr = {{1,2,3}};
            REQUIRE(pr[2] == 3);
            REQUIRE(mod == true);
        }
        SECTION("Atom::charge: wrapped float"){
            bool mod{};
            float f{};
            Atom::PropRef<float> pr = {&f,&mod};
            REQUIRE(pr == 0);
            REQUIRE(mod == false);
            pr = 5.27;
            REQUIRE(f == 5.27f);
            REQUIRE(mod == true);
        }
        //TODO: fixvec should be removed -> property-char
    }

    SECTION( "Vipster::AtomProper: implementation of a free Atom" ) {
        AtomProper ap2{"C",{{1,2,3}}};
        REQUIRE( ap.name == "" );
        REQUIRE( ap.coord == Vec{{0,0,0}} );
        REQUIRE( ap2.name == "C" );
        REQUIRE( ap2.coord == Vec{{1,2,3}} );
        REQUIRE( ap.charge == 0 );
        REQUIRE( ap.charge == ap2.charge );
        REQUIRE( ap.hidden == false );
        REQUIRE( ap.hidden == ap2.hidden );
        ap2.charge = 5.27;
        REQUIRE( ap2.charge == 5.27f );
        AtomProper ap3 = ap2;
        ap2.charge = 6.35;
        REQUIRE( ap2.charge == 6.35f );
        REQUIRE( ap3.charge == 5.27f );
        AtomProper ap4{};
        REQUIRE( ap4.charge == 0 );
        ap4 = ap2;
        ap2.charge = 7.29;
        REQUIRE( ap2.charge == 7.29f );
        REQUIRE( ap4.charge == 6.35f );
    }
//    Atom atp{ap1};
//    Step s{};
//    s.newAtom();
//    Atom ats{s[0]};
}

TEST_CASE( "Vipster::Step", "[step]" ) {
    StepProper s{};

    SECTION( "Step-basics" ) {
        REQUIRE( s.getFmt() == AtomFmt::Bohr );
        REQUIRE( s.getNat() == 0);
        REQUIRE( s.getComment() == "" );
        s.setComment("test");
        REQUIRE( s.getComment() == "test" );
        REQUIRE( s.getCellDim(CdmFmt::Bohr) == 1 );
        REQUIRE( s.getCellDim(CdmFmt::Angstrom) == bohrrad );
        REQUIRE( s.getCellVec() == Mat{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}});
        REQUIRE( s.getNtyp() == 0 );
    }

    auto testStep = [](Step& s){
        AtomProper at1{};
        AtomProper at2{"H", {{1,0,0}}};

        s.newAtom();
        s.newAtom(at1);
        s.newAtom(at2);
        REQUIRE( s[0] == at1 );
        REQUIRE( s[1] == at1 );
        REQUIRE( s[2] != at1 );
        REQUIRE( s[2] == at2 );
        REQUIRE( s.getNat() == 3 );
        REQUIRE( s.getNtyp() == 2 );
        REQUIRE( s.getTypes() == std::set<std::string>{"","H"} );
        s.delAtom(2);
        REQUIRE( s.getNat() == 2 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{""} );
        for(auto& at: s) {
            REQUIRE( at == at1 );
        }
        [at1](const Step& s){
            for(auto& at: s){
                REQUIRE( at == at1 );
            }
        }(s);
        s.newAtoms(3);
        REQUIRE( s.getNat() == 5 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{""} );
        s[0] = at2;
        REQUIRE(s[0] == at2);
    };

    SECTION( "Vipster::StepProper: basic Step class" ) {
        testStep(s);
    }

    SECTION( "Vipster::StepFormatter: formatted interface to StepProper") {
        std::string fmtNames[] = {"asBohr", "asAngstrom", "asCrystal", "asAlat"};
        for(size_t i=0; i<4; ++i){
            SECTION( fmtNames[i] ){
                testStep(s.asFmt((AtomFmt)i));
            }
        }
    }
}
