#include "vipster/step.h"
#include <iostream>

using namespace Vipster;

namespace std {
static ostream& operator<<(ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static ostream& operator<<(ostream& out, const decltype(Step::selection::atom::coord)& v)
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

static ostream& operator<<(ostream& os, const Step::selection::atom& at)
{
    os << "Atom{" << at.name << ", " << at.coord << ", " << at.properties->charge
       << ", " << at.properties->flags << "}";
    return os;
}

}

#include "catch2/catch.hpp"

TEST_CASE("Vipster::detail::Selection", "[select]") {
    Step s{};
    s.setFmt(AtomFmt::Bohr);
    s.newAtom("C");
    s.newAtom("H", Vec{1,0,0});
    s.newAtom("O", Vec{2,2,2});
    s.generateBonds();

    SECTION("by type"){
        auto selC = s.select("type C");
        CHECK(selC.getNat() == 1);
        CHECK(selC[0] == s[0]);
        auto selH = s.select("type H");
        CHECK(selH.getNat() == 1);
        CHECK(selH[0] == s[1]);
        auto selO = s.select("type O");
        CHECK(selO.getNat() == 1);
        CHECK(selO[0] == s[2]);
        auto selCHO = s.select("type [C H O]");
        CHECK(selCHO.getNat() == 3);
        CHECK(selCHO[0] == s[0]);
        CHECK(selCHO[1] == s[1]);
        CHECK(selCHO[2] == s[2]);
    }

    SECTION("by index"){
        auto sel0 = s.select("index 0");
        CHECK(sel0.getNat() == 1);
        CHECK(sel0[0] == s[0]);
        auto sel1 = s.select("index 1");
        CHECK(sel1.getNat() == 1);
        CHECK(sel1[0] == s[1]);
        auto sel2 = s.select("index 2");
        CHECK(sel2.getNat() == 1);
        CHECK(sel2[0] == s[2]);
        auto selRange = s.select("index 0-2");
        CHECK(selRange.getNat() == 3);
        CHECK(selRange[0] == s[0]);
        CHECK(selRange[1] == s[1]);
        CHECK(selRange[2] == s[2]);
        auto selList = s.select("index [0 2]");
        CHECK(selList.getNat() == 2);
        CHECK(selList[0] == s[0]);
        CHECK(selList[1] == s[2]);
        auto selListRange = s.select("index [1-2 0]");
        CHECK(selListRange.getNat() == 3);
        CHECK(selListRange[0] == s[1]);
        CHECK(selListRange[1] == s[2]);
        CHECK(selListRange[2] == s[0]);
        auto selOverdefined = s.select("index [0 0-2]");
        CHECK(selOverdefined.getNat() == 4);
        CHECK(selOverdefined[0] == s[0]);
        CHECK(selOverdefined[1] == s[0]);
        CHECK(selOverdefined[2] == s[1]);
        CHECK(selOverdefined[3] == s[2]);
    }

    SECTION("by pos"){
        s.setCellDim(2, AtomFmt::Bohr);
        s.setCellVec(Mat{Vec{1,1,0}, Vec{0,1,0}, Vec{0,0,1}});
        auto selG = s.select("pos x>0");
        CHECK(selG.getNat() == 2);
        auto selGE = s.select("pos x>= 1");
        CHECK(selGE.getNat() == 2);
        auto selL = s.select("pos y < 2");
        CHECK(selL.getNat() == 2);
        auto selLE = s.select("pos y <= 0");
        CHECK(selLE.getNat() == 2);
    }

    SECTION("by coord"){
        auto selEQ = s.select("coord =1");
        CHECK(selEQ.getNat() == 2);
        auto selGT = s.select("coord >0");
        CHECK(selGT.getNat() == 2);
        auto selLT = s.select("coord <1");
        CHECK(selLT.getNat() == 1);
    }

    SECTION("not"){
        auto notIdxList = s.select("not index [0 2]");
        CHECK(notIdxList.getNat() == 1);
        auto notTypeList = s.select("not type [C O]");
        CHECK(notTypeList.getNat() == 1);
        auto notCoordLT = s.select("not coord <1");
        CHECK(notCoordLT.getNat() == 2);
        auto notPosB = s.select("not pos x>0");
        CHECK(notPosB.getNat() == 1);
    }

    SECTION("complex"){
        auto selOr = s.select("index 0 | type H");
        CHECK(selOr.getNat() == 2);
        auto selNor = s.select("index 0 !| type H");
        CHECK(selNor.getNat() == 1);
        auto selAnd = s.select("index 0-2 & type H");
        CHECK(selAnd.getNat() == 1);
        auto selNand = s.select("index 0-2 !& type H");
        CHECK(selNand.getNat() == 2);
        auto selXor = s.select("not index 2 ^ type [H O]");
        CHECK(selXor.getNat() == 2);
        auto selXnor = s.select("not index 2 !^ type [H O]");
        CHECK(selXnor.getNat() == 1);
        auto ungrouped = s.select("index 0 | index 1 & index 2");
        CHECK(ungrouped.getNat() == 1);
        auto leftgrouped = s.select("(index 0 | index 0-2) & index 2");
        CHECK(leftgrouped.getNat() == 1);
        auto rightgrouped = s.select("index 0 | (index 0-2 & index 2)");
        CHECK(rightgrouped.getNat() == 2);
    }

    SECTION("from selection"){
        auto sel = s.select("index [0 2]");
        CHECK(sel.getNat() == 2);
        auto selSel = sel.select("index 0");
        CHECK(selSel.getNat() == 1);
    }
}
