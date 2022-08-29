#include "vipster/plugins/xyz.h"
#include "vipster/vec.h"
#include <iostream>
#include <sstream>

#include <catch2/catch.hpp>

using namespace Vipster;

TEST_CASE("Vipster::Plugins::XYZ", "[IO][XYZ]") {
    SECTION("ParseEmpty") {
        std::stringstream testInput{""};
        auto [mol, param, data] = Plugins::XYZ.parser("empty", testInput);
        REQUIRE(!param.has_value());
        REQUIRE(data.empty());
        REQUIRE(mol.getNstep() == 0);
        REQUIRE(mol.name == "empty");
    }

    SECTION("ParseSingle") {
        std::stringstream testInput{
            "2\n"
            "comment\n"
            "C 0 0 0\n"
            "C 1.5 0 0\n"
        };
        auto [mol, param, data] = Plugins::XYZ.parser("testname", testInput);
        REQUIRE(!param.has_value());
        REQUIRE(data.empty());
        REQUIRE(mol.getNstep() == 1);
        const auto &s = mol.getStep(0);
        REQUIRE(s.getNat() == 2);
        REQUIRE(s[0].name == "C");
        REQUIRE(s[1].name == "C");
        REQUIRE(s[0].coord == Vec{0,0,0});
        REQUIRE(s[1].coord == Vec{1.5,0,0});
    }
}
