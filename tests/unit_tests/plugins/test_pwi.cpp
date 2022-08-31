#include "vipster/plugins/pwinput.h"
#include "vipster/vec.h"
#include <iostream>
#include <sstream>

#include <catch2/catch.hpp>

using namespace Vipster;
using namespace std::literals::string_literals;
using NameList = std::map<std::string, std::string>;

TEST_CASE("Vipster::Plugins::PWInput", "[IO][PWI]") {
    SECTION("Preset") {
        auto preset = Plugins::PWInput.makePreset();

        // check for default values (index()==1 for NamedEnum)
        REQUIRE(preset.size() == 2);

        REQUIRE(preset.at("atoms").first.index() == 1);
        REQUIRE(std::get<NamedEnum>(preset.at("atoms").first).name() == "Active"s);

        REQUIRE(preset.at("cell").first.index() == 1);
        REQUIRE(std::get<NamedEnum>(preset.at("cell").first).name() == "Active"s);
    }

    SECTION("Parameter") {
        auto param = Plugins::PWInput.makeParam();

        // check for default values (index() == 0: string, == 2: NameList)
        REQUIRE(param.size() == 7);

        const char *namelists[] = {"&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"};
        for (const auto &s: namelists){
            REQUIRE(param.at(s).first.index() == 2);
            REQUIRE(std::get<NameList>(param.at(s).first).empty());
        }

        const char *strings[] = {"PPPrefix", "PPSuffix"};
        for (const auto &s: strings){
            REQUIRE(param.at(s).first.index() == 0);
            REQUIRE(std::get<std::string>(param.at(s).first).empty());
        }
    }

    SECTION("ParseSimple") {
        std::stringstream testInput{
        R"--(
            &CONTROL
             calculation='relax'
             restart_mode='from_scratch'
            /

            &system
             nat=2
             ntyp=1
             celldm(1)=5.0
             ibrav=0
            /

            &electrons
             conv_thr=1.0e-8
            /

            ATOMIC_SPECIES
            C 12.123 C.uspp.pbe.UPF

            ATOMIC_POSITIONS crystal
            C 0.0 0.0 0.0
            C 0.5 0.5 0.5 0 0 1

            K_POINTS gamma

            CELL_PARAMETERS alat
            1 0 0
            0 1 0
            0 0 1
        )--"};

        auto [mol, param, data] = Plugins::PWInput.parser("test", testInput);

        REQUIRE(data.empty());

        // Verify parameter set
        REQUIRE(param.has_value());
        const auto &par = param.value();

        const auto &control = std::get<NameList>(par.at("&CONTROL").first);
        CHECK(control.at("calculation") == "'relax'");
        CHECK(control.at("restart_mode") == "'from_scratch'");

        const auto &system = std::get<NameList>(par.at("&SYSTEM").first);
        // &SYSTEM contains mandatory parameters which are removed from the namelist after parsing
        CHECK(system.find("nat") == system.end());
        CHECK(system.find("ntyp") == system.end());
        CHECK(system.find("celldm(1)") == system.end());
        CHECK(system.find("ibrav") == system.end());

        const auto &electrons = std::get<NameList>(par.at("&ELECTRONS").first);
        CHECK(electrons.at("conv_thr") == "1.0e-8");

        CHECK(std::get<NameList>(par.at("&IONS").first).empty());
        CHECK(std::get<NameList>(par.at("&CELL").first).empty());

        // Verify molecule
        REQUIRE(mol.getNstep() == 1);
        CHECK(mol.name == "test");

        // custom atom properties
        const auto &pte = mol.getPTE();
        CHECK(pte.at("C").PWPP == "C.uspp.pbe.UPF");
        CHECK(pte.at("C").m == 12.123);

        // atoms are in crystal coordinates
        const auto &s = mol.getStep(0);
        REQUIRE(s.hasCell());
        CHECK(s.getCellDim(AtomFmt::Bohr) == 5.0);
        CHECK(s.getCellVec() == Mat{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}});
        CHECK(s[0].name == "C");
        CHECK(s[0].coord == Vec{{0,0,0}});
        CHECK(s[0].properties->flags.none());
        CHECK(s[1].name == "C");
        CHECK(s[1].coord == Vec{{0.5,0.5,0.5}});
        CHECK(s[1].properties->flags[AtomProperties::FixX]);
        CHECK(s[1].properties->flags[AtomProperties::FixY]);
        CHECK(!s[1].properties->flags[AtomProperties::FixZ]);
    }

    SECTION("WriteSimple") {
        std::ostringstream testOutput;

        // configure molecule
        Molecule mol{};
        auto &s = mol.getStep(0);
        s.setFmt(AtomFmt::Alat);
        s.setCellDim(5, AtomFmt::Bohr);
        s.newAtom("C", {0,0,0});
        s.newAtom("C", {0.5,0.5,0.5});
        s[1].properties->flags[AtomProperties::FixZ] = true;

        // configure parameters
        auto param = Plugins::PWInput.makeParam();

        auto &control = std::get<NameList>(param.at("&CONTROL").first);
        control["calculation"] = "'vc-relax'";
        control["nstep"] = "17";
        std::get<std::string>(param.at("PPSuffix").first) = ".uspp.pbe.UPF";

        auto preset = Plugins::PWInput.makePreset();

        bool written = Plugins::PWInput.writer(mol, testOutput, param, preset, 0);

        REQUIRE(written);
        REQUIRE(testOutput.str() ==
R"--(&CONTROL
 calculation = 'vc-relax'
 nstep = 17
/

&SYSTEM
 ibrav = 0
 nat = 2
 ntyp = 1
 celldm(1) = 5
/

&ELECTRONS
/

&IONS
/

&CELL
/

ATOMIC_SPECIES
C    12.01070 C.uspp.pbe.UPF

ATOMIC_POSITIONS alat
C     0.00000000   0.00000000   0.00000000
C     0.50000000   0.50000000   0.50000000 1 1 0

K_POINTS gamma

CELL_PARAMETERS alat
  1.00000000   0.00000000   0.00000000
  0.00000000   1.00000000   0.00000000
  0.00000000   0.00000000   1.00000000
)--");
    }
}
