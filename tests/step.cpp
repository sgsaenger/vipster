#include "catch.hpp"
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

static std::ostream& operator<<(std::ostream& os, Atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.charge
       << ", " << at.fix << ", " << at.hidden << "}";
    return os;
}

TEST_CASE( "Vipster::Atom", "[atom]" ) {

}

TEST_CASE( "Vipster::Step", "[step]" ) {
    StepProper s{};

    REQUIRE( s.getFmt() == AtomFmt::Bohr );
    REQUIRE( s.getNat() == 0);
    REQUIRE( s.getComment() == "" );
    s.setComment("test");
    REQUIRE( s.getComment() == "test" );
    REQUIRE( s.getCellDim() == 1 );
    REQUIRE( s.getCellVec() == Mat{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}});
    REQUIRE( s.getNtyp() == 0 );

    SECTION( "Atom modification through StepProper" ) {
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
        REQUIRE( s.getTypes() == std::set<std::string>{"C","H"} );
        s.delAtom(2);
        REQUIRE( s.getNat() == 2 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{"C"} );
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
        REQUIRE( s.getNtyp() == 2 );
        REQUIRE( s.getTypes() == std::set<std::string>{"C",""} );
    }
}
