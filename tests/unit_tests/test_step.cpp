#include "vipster/step.h"
#include <iostream>

using namespace Vipster;

namespace std {

static ostream& operator<<(ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static ostream& operator<<(ostream& os, const Step::atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.properties->charge
       << ", " << at.properties->flags << "}";
    return os;
}

}

#include "catch2/catch.hpp"

TEST_CASE( "Vipster::Step", "[step]" ) {

    SECTION( "Defaults" ) {
        Step s{};
        // basic fmt-query
        REQUIRE( s.getFmt() == AtomFmt::Angstrom );
        s.setFmt(AtomFmt::Crystal);
        REQUIRE( s.getFmt() == AtomFmt::Crystal );
        // comment
        REQUIRE( s.getComment() == "" );
        s.setComment("test");
        REQUIRE( s.getComment() == "test" );
        // basic cell-functionality
        REQUIRE( s.getCellDim(AtomFmt::Bohr) == Approx(invbohr) );
        REQUIRE( s.getCellDim(AtomFmt::Angstrom) == Approx(1) );
        REQUIRE( s.getCellVec() == Mat{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}});
        // no atom-based properties yet
        REQUIRE( s.getNat() == 0);
        REQUIRE( s.getNtyp() == 0 );
    }

    SECTION( "Basic atom-access" ){
        // empty coordinates don't care for format
        Step s{};
        s.newAtoms(3);
        REQUIRE( s.getNat() == 3 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{""} );
        for(const auto& at: s) {
            REQUIRE( at == s[0] );
            REQUIRE( at.name == "" );
            REQUIRE( at.coord == Vec{} );
            REQUIRE( at.type == s.getPTE().at("") );
            REQUIRE( at.properties == AtomProperties{} );
        }
        s.newAtom("C", {1,2,3}, {1.2, {4,5,6}, AtomProperties::FixX | AtomProperties::Hidden});
        const auto& at = s.back();
        REQUIRE( s.getNat() == 4 );
        REQUIRE( s.getNtyp() == 2 );
        REQUIRE( s.getTypes() == std::set<std::string>{"", "C"} );
        REQUIRE( at.name == "C" );
        REQUIRE( at.coord == Vec{1,2,3} );
        REQUIRE( at.type == s.getPTE().at("C") );
        REQUIRE( at.properties->charge == 1.2 );
        REQUIRE( at.properties->forces == Vec{4,5,6} );
        REQUIRE( at.properties->flags == (AtomProperties::FixX | AtomProperties::Hidden) );
        s.delAtom(3);
        REQUIRE( s.getNat() == 3 );
        REQUIRE( s.getNtyp() == 1 );
        REQUIRE( s.getTypes() == std::set<std::string>{""} );
    }

    SECTION( "Bonds" ){
        Step s{};
        s.newAtom("C", {0,0,0});
        SECTION( "Molecular bonds" ){
            s.newAtom("H", {1,0,0});
            s.generateBonds();
            REQUIRE( s.getBonds().size() == 1 );
            const auto &bond = s.getBonds()[0];
            REQUIRE( bond.at1 == 0 );
            REQUIRE( bond.at2 == 1 );
            REQUIRE( bond.diff == DiffVec{0,0,0} );
            REQUIRE( bond.dist == invbohr );
            REQUIRE( bond.type == nullptr );
        }
        SECTION( "Periodic bonds" ){
            s.setFmt(AtomFmt::Crystal);
            s.newAtom("C", {0.5,0.5,0.5});
            s.setCellDim(3, AtomFmt::Bohr, true);
            s.generateBonds();
            REQUIRE( s.getBonds().size() == 8 );
            for(const auto &bond: s.getBonds()){
                REQUIRE( bond.at1 == 0 );
                REQUIRE( bond.at2 == 1 );
                REQUIRE( bond.dist == std::sqrt(3*1.5*1.5) );
                REQUIRE( bond.type == nullptr );
            }
            REQUIRE( s.getBonds()[0].diff == DiffVec{ 0, 0, 0} );
            REQUIRE( s.getBonds()[1].diff == DiffVec{-1, 0, 0} );
            REQUIRE( s.getBonds()[2].diff == DiffVec{ 0,-1, 0} );
            REQUIRE( s.getBonds()[3].diff == DiffVec{-1,-1, 0} );
            REQUIRE( s.getBonds()[4].diff == DiffVec{ 0, 0,-1} );
            REQUIRE( s.getBonds()[5].diff == DiffVec{-1, 0,-1} );
            REQUIRE( s.getBonds()[6].diff == DiffVec{ 0,-1,-1} );
            REQUIRE( s.getBonds()[7].diff == DiffVec{-1,-1,-1} );
        }
        SECTION( "Manual bonds" ){
            s.newAtom("C", {1,2,3});
            REQUIRE( s.getBonds().size() == 0 );
            s.addBond( 0, 1, {1,2,3}, "CC" );
            REQUIRE( s.getBonds().size() == 1 );
            REQUIRE( s.getBonds()[0].type );
            REQUIRE( s.getBonds()[0].type->first == "CC" );
            s.setBondType(0, "");
            REQUIRE( s.getBonds()[0].type == nullptr );
            s.delBond(0);
            REQUIRE( s.getBonds().size() == 0 );
        }
    }
}
