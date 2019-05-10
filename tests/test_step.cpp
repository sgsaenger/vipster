#include "step.h"
#include "configfile.h"
#include <iostream>

using namespace Vipster;

namespace std {

static ostream& operator<<(ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static ostream& operator<<(ostream& os, const Atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.properties->charge
       << ", " << at.properties->flags << "}";
    return os;
}

}

#include "catch.hpp"

TEST_CASE( "Vipster::Step", "[step]" ) {
    Mat cv = {{{{1,2,3}},{{0,1,0}},{{0,0,1.5}}}};
    std::string fmtNames[nAtFmt] = {"Bohr", "Angstrom", "Crystal", "Alat"};
    // s will be checked

    SECTION( "Defaults, Format, Cell" ) {
        Step s{};
        // basic fmt-query
        REQUIRE( s.getFmt() == AtomFmt::Bohr );
        s.setFmt(AtomFmt::Crystal);
        REQUIRE( s.getFmt() == AtomFmt::Crystal );
        // comment
        REQUIRE( s.getComment() == "" );
        s.setComment("test");
        REQUIRE( s.getComment() == "test" );
        // basic cell-functionality
        REQUIRE( s.getCellDim(CdmFmt::Bohr) == Approx(1) );
        REQUIRE( s.getCellDim(CdmFmt::Angstrom) == Approx(bohrrad) );
        s.setCellDim(1, CdmFmt::Angstrom);
        REQUIRE( s.getCellDim(CdmFmt::Bohr) == Approx(invbohr) );
        REQUIRE( s.getCellDim(CdmFmt::Angstrom) == Approx(1) );
        REQUIRE( s.getCellVec() == Mat{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}});
        s.setCellVec(cv);
        REQUIRE( s.getCellVec() == cv );
        // no atom-based properties yet
        REQUIRE( s.getNat() == 0);
        REQUIRE( s.getNtyp() == 0 );
        REQUIRE( s.getNbond() == 0 );
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
        }
    }

    SECTION( "Formatted coordinates") {
        Step s{};
        s.setCellDim(5, CdmFmt::Bohr);
        s.setCellVec(cv);
        s.newAtom();
        s.newAtom("H", {{1,2,3}});
        // s_comp will be checked against
        Step s_comp[nAtFmt];
        for(size_t i=0; i<nAtFmt; ++i){
            s_comp[i].setFmt(static_cast<AtomFmt>(i));
            s_comp[i].setCellDim(5, CdmFmt::Bohr);
            s_comp[i].setCellVec(cv);
            // generic comparison targets
            s_comp[i].newAtom();
            s_comp[i].newAtom("H", {{1,2,3}});
        }
        // Bohr:
        s_comp[0].newAtom("H", {{1,2,3}});
        s_comp[0].newAtom("H", {{1*bohrrad, 2*bohrrad, 3*bohrrad}});
        s_comp[0].newAtom("H", {{1/5.f, 0, 0}});
        s_comp[0].newAtom("H", {{1/5.f, 2/5.f, 3/5.f}});
        // Angstrom:
        s_comp[1].newAtom("H", {{1*invbohr, 2*invbohr, 3*invbohr}});
        s_comp[1].newAtom("H", {{1, 2, 3}});
        s_comp[1].newAtom("H", {{1*invbohr/5.f, 0, 0}});
        s_comp[1].newAtom("H", {{1*invbohr/5.f, 2*invbohr/5.f, 3*invbohr/5.f}});
        // Crystal:
        s_comp[2].newAtom("H", {{5, 20, 37.5f}});
        s_comp[2].newAtom("H", {{5*bohrrad, 20*bohrrad, 37.5f*bohrrad}});
        s_comp[2].newAtom("H", {{1, 2, 3}});
        s_comp[2].newAtom("H", {{1, 4, 7.5f}});
        // Alat:
        s_comp[3].newAtom("H", {{5, 10, 15}});
        s_comp[3].newAtom("H", {{5*bohrrad, 10*bohrrad, 15*bohrrad}});
        s_comp[3].newAtom("H", {{1, 0, 0}});
        s_comp[3].newAtom("H", {{1, 2, 3}});

        for(size_t i=0; i<nAtFmt; ++i){
        SECTION( fmtNames[i] ){
            auto fmt_i = static_cast<AtomFmt>(i);
            s.modScale(fmt_i);

            REQUIRE( s[0] == s_comp[i][0] );
            REQUIRE( s[1] == s_comp[i][1] );
            REQUIRE( s[0] != s_comp[i][1] );
            REQUIRE( s[1] != s_comp[i][0] );
            REQUIRE( s.getNat() == 2 );
            REQUIRE( s.getNtyp() == 2 );
            REQUIRE( s.getTypes() == std::set<std::string>{"","H"} );
            s[0] = s[1];
            for (auto& at: s) {
                REQUIRE( at == s[0]);
            }
            for (size_t j=0; j<nAtFmt; ++j){
                auto fmt_j = static_cast<AtomFmt>(j);
                auto f = s.asFmt(fmt_j);
                f.evaluateCache();
                if(fmt_j != fmt_i){
                    REQUIRE( f[1] != s[1] );
                } else {
                    REQUIRE( f[1] == s[1] );
                }
                CHECK( f[1] == s_comp[i][2+j]);
            }
            s.delAtom(1);
            s.delAtom(0);
            REQUIRE( s.getNat() == 0 );
            REQUIRE( s.getNtyp() == 0 );
            REQUIRE( s.getTypes() == std::set<std::string>{} );
        }
        }
    }

    SECTION( "Cell modification" ){
        Vec v = {{1,2,3}};
        Step s{};
        s.setCellDim(5, CdmFmt::Bohr);
        s.setCellVec(cv);
        s.newAtom("H", v);
        for(size_t i=0; i<nAtFmt; ++i){
        SECTION( fmtNames[i] ){
            auto fmt_i = static_cast<AtomFmt>(i);
            s.modScale(fmt_i);
            SECTION("CellDim"){
                Step sd{}; // alternating unscaled and scaled for each format
                sd.setCellDim(5, CdmFmt::Bohr);
                sd.setCellVec(cv);
                // Bohr:
                sd.newAtom("H", v);
                sd.newAtom("H", v*1.2f);
                // Angstrom:
                sd.newAtom("H", v*bohrrad);
                sd.newAtom("H", v*1.2f*bohrrad);
                // Crystal:
                sd.newAtom("H", Vec{{1,0,0}}/5);
                sd.newAtom("H", Vec{{1,0,0}}/6);
                // Alat:
                sd.newAtom("H", v/5);
                sd.newAtom("H", v/6);
                // scale it so `v` is at the position of proper fmt
                sd.modScale(fmt_i);
                sd.setFmt(AtomFmt::Bohr);
                for(size_t j=0; j<nAtFmt; ++j){
                    INFO(fmtNames[j])
                    auto fmt_j = static_cast<AtomFmt>(j);
                    auto f = s.asFmt(fmt_j);
                    f.evaluateCache();
                    bool rel_i = fmt_i>=AtomFmt::Crystal;
                    bool abs_i = !rel_i;
                    bool rel_j = fmt_j>=AtomFmt::Crystal;
                    bool abs_j = !rel_j;
                    CHECK(s[0] == sd[2*i]);
                    CHECK(f[0] == sd[2*j]);
                    f.setCellDim(6, CdmFmt::Bohr, false);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sd[2*i+rel_i]);
                    CHECK(f[0] == sd[2*j+rel_j]);
                    f.setCellDim(5, CdmFmt::Bohr, false);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sd[2*i]);
                    CHECK(f[0] == sd[2*j]);
                    f.setCellDim(6, CdmFmt::Bohr, true);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sd[2*i+abs_i]);
                    CHECK(f[0] == sd[2*j+abs_j]);
                    f.setCellDim(5, CdmFmt::Bohr, true);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sd[2*i]);
                    CHECK(f[0] == sd[2*j]);
                }
            }
            SECTION("CellVec"){
                Step sv{}; // alternating unscaled and scaled for each format
                bool rel_i = fmt_i == AtomFmt::Crystal;
                bool abs_i = !rel_i;
                Mat cv2 = {{{{1,0,4}},{{0,1,0}},{{0,0,2}}}};
                Vec v2 = {{1,0,4}};
                sv.setCellDim(5, CdmFmt::Bohr);
                sv.setCellVec(cv);
                if(abs_i){
                    sv.setFmt(fmt_i);
                    // Bohr:
                    sv.newAtom("H", v);
                    sv.newAtom("H", v2);
                    // Angstrom:
                    sv.newAtom("H", v*bohrrad);
                    sv.newAtom("H", v2*bohrrad);
                    // Crystal:
                    sv.newAtom("H", Vec{{1,0,0}}/5);
                    sv.newAtom("H", Vec{{1,2,-0.5}}/5);
                    // Alat:
                    sv.newAtom("H", v/5);
                    sv.newAtom("H", v2/5);
                    // scale it so `v` is at the position of proper fmt
                    sv.setFmt(AtomFmt::Bohr);
                } else {
                    // Bohr:
                    sv.newAtom("H", {{5,20,37.5}});
                    sv.newAtom("H", {{5,10,50}});
                    // Angstrom:
                    sv.newAtom("H", Vec{{5,20,37.5}}*bohrrad);
                    sv.newAtom("H", Vec{{5,10,50}}*bohrrad);
                    // Crystal:
                    sv.newAtom("H", {{1,2,3}});
                    sv.newAtom("H", {{1,4,1.75}});
                    // Alat:
                    sv.newAtom("H", {{1,4,7.5}});
                    sv.newAtom("H", {{1,2,10}});
                }
                for(size_t j=0; j<nAtFmt; ++j){
                    INFO(fmtNames[j])
                    auto fmt_j = static_cast<AtomFmt>(j);
                    auto f = s.asFmt(fmt_j);
                    f.evaluateCache();
                    bool rel_j = fmt_j == AtomFmt::Crystal;
                    bool abs_j = !rel_j;
                    CHECK(s[0] == sv[2*i]);
                    CHECK(f[0] == sv[2*j]);
                    f.setCellVec(cv2, false);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sv[2*i+rel_i]);
                    CHECK(f[0] == sv[2*j+rel_j]);
                    f.setCellVec(cv, false);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sv[2*i]);
                    CHECK(f[0] == sv[2*j]);
                    f.setCellVec(cv2, true);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sv[2*i+abs_i]);
                    CHECK(f[0] == sv[2*j+abs_j]);
                    f.setCellVec(cv, true);
                    f.evaluateCache();
                    s.evaluateCache();
                    CHECK(s[0] == sv[2*i]);
                    CHECK(f[0] == sv[2*j]);
                }
            }
        }
        }
    }
}
