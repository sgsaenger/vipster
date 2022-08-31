#include "vipster/plugins/xyz.h"
#include "vipster/vec.h"
#include <iostream>
#include <sstream>

#include <catch2/catch.hpp>

using namespace Vipster;
using namespace std::literals::string_literals;

TEST_CASE("Vipster::Plugins::XYZ", "[IO][XYZ]") {
    SECTION("Preset") {
        auto preset = Plugins::XYZ.makePreset();

        // check for default values (index()==1 for NamedEnum)
        REQUIRE(preset.size() == 2);

        REQUIRE(preset.find("filemode") != preset.end());
        REQUIRE(preset.at("filemode").first.index() == 1);
        REQUIRE(std::get<NamedEnum>(preset.at("filemode").first).name() == "Step"s);

        REQUIRE(preset.find("atomdata") != preset.end());
        REQUIRE(preset.at("atomdata").first.index() == 1);
        REQUIRE(std::get<NamedEnum>(preset.at("atomdata").first).name() == "None"s);
    }

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
        REQUIRE(mol.name == "testname");

        const auto &s = mol.getStep(0);
        REQUIRE(s.getNat() == 2);
        REQUIRE(s.getComment() == "comment");
        REQUIRE(s[0].name == "C");
        REQUIRE(s[1].name == "C");
        REQUIRE(s[0].coord == Vec{0,0,0});
        REQUIRE(s[1].coord == Vec{1.5,0,0});
        REQUIRE(!s.hasCell());
    }

    SECTION("ParseTrajec") {
        std::stringstream testInput{
            "2\ncomment\nC 0 0 0\nC 1.5 0 0\n"
            "2\ncomment\nC 0 0 0\nC 1.5 0 0\n"
            "\n"
            "2\ncomment\nC 0 0 0\nC 1.5 0 0\n"
        };

        auto [mol, param, data] = Plugins::XYZ.parser("testname", testInput);
        REQUIRE(!param.has_value());
        REQUIRE(data.empty());

        REQUIRE(mol.getNstep() == 3);
        REQUIRE(mol.name == "testname");

        for(const auto &s: mol.getSteps()){
            REQUIRE(s.getNat() == 2);
            REQUIRE(s.getComment() == "comment");
            REQUIRE(s[0].name == "C");
            REQUIRE(s[1].name == "C");
            REQUIRE(s[0].coord == Vec{0,0,0});
            REQUIRE(s[1].coord == Vec{1.5,0,0});
            REQUIRE(!s.hasCell());
        }
    }

    SECTION("ParseCell") {
        std::stringstream testInput{
            "2\n"
            "comment\n"
            "C 0 0 0\n"
            "C 1.5 0 0\n"
            "\n"
            "3 0 0\n"
            "0 3 0\n"
            "0 0 3\n"
        };

        auto [mol, param, data] = Plugins::XYZ.parser("testname", testInput);
        REQUIRE(!param.has_value());
        REQUIRE(data.empty());

        REQUIRE(mol.getNstep() == 1);
        REQUIRE(mol.name == "testname");

        const auto &s = mol.getStep(0);
        REQUIRE(s.getNat() == 2);
        REQUIRE(s.getComment() == "comment");
        REQUIRE(s[0].name == "C");
        REQUIRE(s[1].name == "C");
        REQUIRE(s[0].coord == Vec{0,0,0});
        REQUIRE(s[1].coord == Vec{1.5,0,0});
        REQUIRE(s.hasCell());
        REQUIRE(s.getCellVec() == Mat{{{{3,0,0}}, {{0,3,0}}, {{0,0,3}}}});
        REQUIRE(s.getCellDim(AtomFmt::Angstrom) == 1);
    }

    SECTION("WriteSingle") {
        std::ostringstream testOutput;
        Molecule mol{};
        auto &s = mol.getStep(0);
        s.setComment("comment");
        s.newAtom("C", {0,0,0});
        s.newAtom("C", {1.5,0,0});

        bool written = Plugins::XYZ.writer(mol, testOutput, {}, Plugins::XYZ.makePreset(), 0);

        REQUIRE(written);
        REQUIRE(testOutput.str() ==
            "2\n"
            "comment\n"
            "C      0.00000    0.00000    0.00000\n"
            "C      1.50000    0.00000    0.00000\n"
        );
    }

    SECTION("WriteTrajec") {
        std::ostringstream testOutput;
        Molecule mol{"trajec", 3};
        for(auto &s: mol.getSteps()){
            s.setComment("comment");
            s.newAtom("C", {0,0,0});
            s.newAtom("C", {1.5,0,0});
        }

        auto preset = Plugins::XYZ.makePreset();
        std::get<NamedEnum>(preset.at("filemode").first) = "Trajec";

        bool written = Plugins::XYZ.writer(mol, testOutput, {}, preset, 0);

        REQUIRE(written);
        REQUIRE(testOutput.str() ==
            "2\ncomment\nC      0.00000    0.00000    0.00000\nC      1.50000    0.00000    0.00000\n"
            "2\ncomment\nC      0.00000    0.00000    0.00000\nC      1.50000    0.00000    0.00000\n"
            "2\ncomment\nC      0.00000    0.00000    0.00000\nC      1.50000    0.00000    0.00000\n"
        );
    }

    SECTION("WriteCell") {
        std::ostringstream testOutput;
        Molecule mol{};
        auto &s = mol.getStep(0);
        s.setComment("comment");
        s.newAtom("C", {0,0,0});
        s.newAtom("C", {1.5,0,0});
        s.setCellVec(Mat{{{{3,0,0}}, {{0,3,0}}, {{0,0,3}}}});

        auto preset = Plugins::XYZ.makePreset();
        std::get<NamedEnum>(preset.at("filemode").first) = "Cell";

        bool written = Plugins::XYZ.writer(mol, testOutput, {}, preset, 0);

        REQUIRE(written);
        REQUIRE(testOutput.str() ==
            "2\n"
            "comment\n"
            "C      0.00000    0.00000    0.00000\n"
            "C      1.50000    0.00000    0.00000\n"
            "\n"
            "3.00000 0.00000 0.00000\n"
            "0.00000 3.00000 0.00000\n"
            "0.00000 0.00000 3.00000\n"
        );
    }
}
