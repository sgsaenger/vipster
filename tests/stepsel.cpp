#include "stepsel.h"
#include "catch.hpp"
using namespace Vipster;

TEST_CASE("Vipster::StepSel", "[step]") {
    StepProper s{};
    s.setFmt(AtomFmt::Bohr);
    s.newAtom("C");
    s.newAtom("H", Vec{1,0,0});
    s.newAtom("O", Vec{2,2,2});
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
        CHECK(selListRange[0] == s[0]);
        CHECK(selListRange[1] == s[1]);
        CHECK(selListRange[2] == s[2]);
        auto selOverdefined = s.select("index [0 0-2]");
        CHECK(selOverdefined.getNat() == 3);
        CHECK(selOverdefined[0] == s[0]);
        CHECK(selOverdefined[1] == s[1]);
        CHECK(selOverdefined[2] == s[2]);
    }
    SECTION("by pos"){
        s.setCellDim(2, CdmFmt::Bohr);
        s.setCellVec(Mat{Vec{1,1,0}, Vec{0,1,0}, Vec{0,0,1}});
        auto selA = s.select("pos xa>0");
        CHECK(selA.getNat() == 2);
        auto selB = s.select("pos xb>0");
        CHECK(selB.getNat() == 2);
        auto selC = s.select("pos y c < 0");
        CHECK(selC.getNat() == 1);
        auto selD = s.select("pos xd > 0.5");
        CHECK(selD.getNat() == 1);
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
        auto notPosB = s.select("not pos xb>0");
        CHECK(notPosB.getNat() == 1);
    }
}
