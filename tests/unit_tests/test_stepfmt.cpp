#include "vipster/step.h"
#include <iostream>

using namespace Vipster;

namespace std {

static ostream& operator<<(ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static ostream& operator<<(ostream& out, const decltype(Step::atom::coord)& v)
{
    out << static_cast<const Vec&>(v);
    return out;
}

static ostream& operator<<(ostream& out, const decltype(Step::formatter::atom::coord)& v)
{
    out << static_cast<const Vec&>(v);
    return out;
}

static ostream& operator<<(ostream& os, const Step::const_atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.properties->charge
       << ", " << at.properties->flags << "}";
    return os;
}

static ostream& operator<<(ostream& os, const Step::formatter::atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.properties->charge
       << ", " << at.properties->flags << "}";
    return os;
}

}

#include "catch2/catch.hpp"

TEST_CASE( "Vipster::detail::Formatter", "[format]"){
    Mat cv = {{{{1,2,3}},{{0,1,0}},{{0,0,1.5}}}};
    std::string fmtNames[4] = {"Crystal", "Alat", "Angstrom"};

    SECTION( "Formatted coordinates") {
        Step s{};
        s.setCellDim(5, AtomFmt::Angstrom);
        s.setCellVec(cv);
        // add formatted atoms
        s.setFmt(AtomFmt::Crystal);
        s.newAtom("H", {{1,2,3}});
        s.setFmt(AtomFmt::Alat);
        s.newAtom("H", {{1,2,3}});
        s.setFmt(AtomFmt::Angstrom);
        s.newAtom("H", {{1,2,3}});
        // reset s to Angstrom
        s.setFmt(AtomFmt::Angstrom);
        for(size_t i=0; i<3; ++i){
        SECTION("comparison "+fmtNames[i]){
            auto fmt = static_cast<AtomFmt>(i-2);
            Vec vec_comp[3][3] = {{{1,2,3}, {1,0,0}, {0.2,0,0}}, // crystal
                                  {{1,4,7.5}, {1,2,3}, {0.2,0.4,0.6}}, // alat
                                  {{5, 20, 37.5}, {5,10,15}, {1,2,3}}}; // angstrom
            for(size_t j=0; j<3; ++j){
                // atoms are unchanged upon insertion/assignment -> coordinates (possibly) modified
                REQUIRE( s.asFmt(fmt)[j] == s[j] );
                // coords are assigned as-is to the formatter -> resulting position depends on format
                REQUIRE( s.asFmt(fmt)[j].coord == vec_comp[i][j] );
            }
        }
        }
        for(size_t i=0; i<3; ++i){
        SECTION("coord assignment: "+fmtNames[i]){
            auto fmt = static_cast<AtomFmt>(i-2);
            s.newAtom(s[i]);
            for(size_t j=0; j<3; ++j){
                REQUIRE( s.back().coord == s[i].coord );
                REQUIRE( s.asFmt(fmt).back().coord == Vec{1,2,3} );
            }
        }
        }
        for(size_t i=0; i<3; ++i){
        SECTION("atom assignment: "+fmtNames[i]){
            auto fmt = static_cast<AtomFmt>(i-2);
            s.newAtom();
            for(size_t j=0; j<3; ++j){
                s.asFmt(fmt).back() = s[i];
                REQUIRE( s.back().coord == s[i].coord );
            }
            for(size_t j=0; j<3; ++j){
                s.back() = s.asFmt(fmt)[i];
                REQUIRE( s.back().coord == s[i].coord );
            }
        }
        }
    }

}
