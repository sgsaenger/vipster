#include "json.h"

#include <nlohmann/json.hpp>
#include <fmt/format.h>

using namespace Vipster;
using namespace nlohmann;

const std::map<std::string, int> str2fmt = {{"crystal",-2}, {"alat",-1}, {"angstrom",0}, {"bohr",1}};

static Preset makePreset()
{
    return {&Plugins::JSON, {
            {"trajectory", {false, "Write complete trajectory instead of single step"}},
            {"bonds", {false, "Export bonds"}},
            {"cell", {false, "Write with cell"}},
            {"fmt", {NamedEnum{4, {"crystal", "alat", "angstrom", "bohr", "active"}},
                "Active: use the current Step's active atom format\n"
                "Else: enforce the selected format"}},
            {"indent", {false, "Pretty-print output so not everything is on one line"}}
           }};
}

IOTuple JSONParser(const std::string& name, std::istream &file)
{
    Molecule m{name, 0};
    json jj{};
    try{
        file >> jj;
        if(!jj.is_array()) throw IOError("JSON-Parse: not a list of steps");
        for(const auto& j: jj){
            if (auto atoms = j.find("atoms"); atoms != j.end()) {
                auto &s = m.newStep();
                // handle format
                auto fmt = AtomFmt::Angstrom;
                if (auto f = j.find("fmt"); f != j.end()) {
                    fmt = static_cast<AtomFmt>(str2fmt.at(f->get<std::string>()));
                }
                s.setFmt(fmt);
                // parse atoms
                if(!atoms->is_array()){
                    throw IOError(fmt::format("JSON-Parser: Not a valid atom-array: {}", atoms->dump()));
                }
                for (const auto &at: *atoms) {
                    auto name = at.find("name");
                    if(name == at.end() || !name->is_string()){
                        throw IOError(fmt::format("JSON-Parser: Invalid atom (missing or invalid \"name\"): {}", at.dump()));
                    }
                    auto coord = at.find("coord");
                    if(coord == at.end() || !coord->is_array() || !(coord->size() == 3)){
                        throw IOError(fmt::format("JSON-Parser: Invalid atom (missing or invalid \"coord\"): {}", at));
                    }
                    s.newAtom(name->get<std::string>(), {(*coord)[0].get<double>(), (*coord)[1].get<double>(), (*coord)[2].get<double>()});
                }
                // parse bonds
                if(auto bonds = j.find("bonds"); bonds != j.end() && bonds->is_array()){
                    for(const auto &bond: *bonds){
                        s.addBond(bond[0], bond[1]);
                    }
                }
                // parse cell
                if(auto cell = j.find("cell"); cell != j.end()){
                    auto c = *cell;
                    s.enableCell(true);
                    s.setCellDim(c["dimension"].get<double>(), s.getFmt() == AtomFmt::Bohr ? AtomFmt::Bohr : AtomFmt::Angstrom);
                    s.setCellVec(c["vectors"].get<Mat>());
                }
            }
        }
    }catch(nlohmann::detail::exception& e){
        throw Error(e.what());
    }
    return {std::move(m), std::optional<Parameter>{}, DataList{}};
}

bool JSONWriter(const Molecule &m, std::ostream &file,
                const std::optional<Parameter> &,
                const std::optional<Preset> &p,
                size_t index)
{
    if(!p || p->getFmt() != &Plugins::JSON){
        throw IOError("JSON-Writer: needs suitable IO preset");
    }
    auto pp = *p;
    // setup json
    auto jj = json::array();
    // iterate over relevant steps
    auto trajec = std::get<bool>(pp.at("trajectory").first);
    auto beg = trajec ? m.getSteps().begin() : std::next(m.getSteps().begin(), index);
    auto end = trajec ? m.getSteps().end() : std::next(beg, 1);
    for (auto it = beg; it != end; ++it) {
        jj.push_back(json::object());
        auto &j = jj.back();
        // get step in correct format
        auto fmt = std::get<NamedEnum>(pp.at("fmt").first);
        if (fmt.value() == 4) fmt = static_cast<int>(it->getFmt())+2;
        j["fmt"] = fmt.name();
        const auto& step = it->asFmt(static_cast<AtomFmt>(fmt.value()-2));
        // store atoms
        auto &atoms = j["atoms"] = json::array();
        for (const auto &at: step) {
            atoms.push_back({
                {"name", at.name.c_str()},
                {"coord", (const Vec&)at.coord}
            });
        }
        // store bonds
        if(std::get<bool>(pp.at("bonds").first)){
            auto& bonds = j["bonds"] = json::array();
            for (const auto& b: step.getBonds()) {
                bonds.push_back({b.at1, b.at2});
            }
        }
        // store cell
        if(std::get<bool>(pp.at("cell").first) && step.hasCell()){
            auto& cell = j["cell"] = json::object();
            cell["dimension"] = step.getCellDim(step.getFmt() == AtomFmt::Bohr ?
                                                AtomFmt::Bohr :
                                                AtomFmt::Angstrom);
            cell["vectors"] = step.getCellVec();
        }
    }
    file << jj.dump(std::get<bool>(pp.at("indent").first) ? 0 : -1);
    return true;
}

const Plugin Plugins::JSON =
{
    "json",
    "json",
    "json",
    &JSONParser,
    &JSONWriter,
    nullptr,
    &makePreset
};
