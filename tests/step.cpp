#include "step.h"
#include <iostream>

using namespace Vipster;

namespace std {

static std::ostream& operator<<(std::ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static std::ostream& operator<<(std::ostream& os, const Atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.charge
       << ", " << at.properties << "}";
    return os;
}

}

#include "catch.hpp"

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
        s.newAtom();
//        s.newAtom(at1);
//        s.newAtom(at2);
//        REQUIRE( s[0] == at1 );
//        REQUIRE( s[1] == at1 );
//        REQUIRE( s[2] != at1 );
//        REQUIRE( s[2] == at2 );
        REQUIRE( s.getNat() == 3 );
        REQUIRE( s.getNtyp() == 2 );
        REQUIRE( s.getTypes() == std::set<std::string>{"","H"} );
        s.delAtom(2);
        REQUIRE( s.getNat() == 2 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{""} );
        for(auto& at: s) {
//            REQUIRE( at == at1 );
        }
//        [at1](const Step& s){
//            for(auto& at: s){
//                REQUIRE( at == at1 );
//            }
//        }(s);
        s.newAtoms(3);
        REQUIRE( s.getNat() == 5 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{""} );
//        s[0] = at2;
//        REQUIRE(s[0] == at2);
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
