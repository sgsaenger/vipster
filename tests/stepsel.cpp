#include "stepsel.h"
#include "catch.hpp"
using namespace Vipster;

TEST_CASE("Vipster::StepSel", "[step]") {
    StepProper s{};
    s.newAtom("C");
    s.newAtom("H");
    s.newAtom("O");
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
