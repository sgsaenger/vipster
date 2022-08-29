#include "json.h"

#include <nlohmann/json.hpp>
#include <fmt/format.h>

using namespace Vipster;
using namespace nlohmann;

const std::map<std::string, AtomFmt> str2fmt = {{"crystal", AtomFmt::Crystal},
                                                {"alat", AtomFmt::Alat},
                                                {"angstrom", AtomFmt::Angstrom},
                                                {"bohr", AtomFmt::Bohr}};

static Preset makePreset()
{
    return {&Plugins::JSON, {
            {"trajectory", {false, "Write complete trajectory instead of single step"}},
            {"bonds", {false, "Export bonds"}},
            {"cell", {false, "Write with cell"}},
            {"fmt", {NamedEnum{0, {"active", "crystal", "alat", "angstrom", "bohr"}},
                "Active: use the current Step's active atom format\n"
                "Else: enforce the selected format"}},
            {"indent", {false, "Pretty-print output so not everything is on one line"}},
            {"elements", {NamedEnum{1, {"none", "custom", "used", "all"}},
                "Store atom types\n"
                "Custom: save only types not in the global table\n"
                "Used: save all types that are used by any atom\n"
                "All: save all types in the molecules' own table"}}
           }};
}

// forward declaration as the json-API is private
namespace Vipster{
    void to_json(json& j, const Element& p);
    void from_json(const json& j, Element& p);
}

void parseStep(const json &j, Molecule &m)
{
    if (auto atoms = j.find("atoms"); atoms != j.end()) {
        auto &s = m.newStep();
        // handle format
        auto fmt = AtomFmt::Angstrom;
        if (auto f = j.find("fmt"); f != j.end()) {
            auto ff = f->get<std::string>();
            auto pos = str2fmt.find(ff);
            if(pos == str2fmt.end()){
                throw IOError(fmt::format("JSON-Parser: Invalid atom format: {} (Valid: [crystal, alat, angstrom, bohr]", ff));
            }
            fmt = pos->second;
        }
        s.setFmt(fmt, false);
        // parse cell
        if(auto cell = j.find("cell"); cell != j.end()){
            auto c = *cell;
            s.enableCell(true);
            s.setCellDim(c["dimension"].get<double>(), s.getFmt() == AtomFmt::Bohr ? AtomFmt::Bohr : AtomFmt::Angstrom);
            s.setCellVec(c["vectors"].get<Mat>());
        }
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
                // check if atoms are valid
                auto at_pos = bond.find("atoms");
                if(at_pos == bond.end()){
                    throw IOError(fmt::format("JSON-Parser: Missing atom ids in bond {}", bond));
                }
                const auto &atoms = at_pos.value();
                if(!(atoms.is_array() && (atoms.size() == 2) &&
                     atoms[0].is_number_unsigned() && atoms[1].is_number_unsigned() &&
                     (atoms[0].get<size_t>() < s.getNat()) && (atoms[1].get<size_t>() < s.getNat()))){
                    throw IOError(fmt::format("JSON-Parser: Invalid atom ids in bond {}", bond));
                }
                // parse optional pbc information
                DiffVec diff{};
                auto pbc_pos = bond.find("pbc");
                if(pbc_pos != bond.end()){
                    const auto &pbc = pbc_pos.value();
                    if(!(pbc.is_array() && pbc.size() == 3 &&
                         pbc[0].is_number() && pbc[1].is_number() && pbc[2].is_number())){
                        throw IOError(fmt::format("JSON-Parser: Invalid periodic bond {}", bond));
                    }
                    diff = {pbc[0].get<int>(), pbc[1].get<int>(), pbc[2].get<int>()};
                }
                // parse optional type string
                std::string type{};
                auto t = bond.find("type");
                if(t != bond.end()){
                    type = t.value().get<std::string>();
                }
                // register bond
                s.addBond(atoms[0], atoms[1], diff, type);
            }
        }
    }
}

IOTuple JSONParser(const std::string& name, std::istream &file)
{
    Molecule m{name, 0};
    json jj{};
    try{
        file >> jj;
        if(!jj.is_object()) throw IOError("JSON-Parser: not a valid json object");
        auto pte = jj.find("elements");
        if(pte != jj.end() && pte->is_object()){
            auto &p = m.getPTE();
            for(const auto &[key, val]: pte->items()){
                p[key] = val;
            }
        }
        auto mol = jj.find("molecule");
        if(mol == jj.end()) throw IOError("JSON-Parser: no valid molecule data");
        if(mol->is_object()){
            parseStep(mol.value(), m);
        }else if(mol->is_array()){
            for(const auto& j: *mol){
                parseStep(j, m);
            }
        }else{
            throw IOError(fmt::format("JSON-Parser: invalid molecule data {}", mol.value()));
        }
    }catch(nlohmann::detail::exception& e){
        throw IOError(e.what());
    }
    return {std::move(m), std::optional<Parameter>{}, DataList{}};
}

void writeStep(json &j, const Step &s, NamedEnum fmt, bool write_cell, bool write_bonds){
    // get step in correct format
    if(fmt.value() == 0) fmt = static_cast<int>(s.getFmt()+3);
    j["fmt"] = fmt.name();
    const auto& step = s.asFmt(static_cast<AtomFmt>(fmt.value()-3));
    // store atoms
    auto &atoms = j["atoms"] = json::array();
    for (const auto &at: step) {
        atoms.push_back({
            {"name", at.name.c_str()},
            {"coord", static_cast<const Vec&>(at.coord)}
        });
    }
    // store bonds
    if(write_bonds){
        auto& bonds = j["bonds"] = json::array();
        DiffVec ref0{0,0,0};
        for (const auto& b: step.getBonds()) {
            if(b.diff != ref0){
                if(write_cell){
                    bonds.push_back({{"atoms", {b.at1, b.at2}},
                                     {"pbc", b.diff}});
                }
            }else{
                bonds.push_back({{"atoms", {b.at1, b.at2}}});
            }
            if(b.type){
                bonds.back()["type"] = b.type->first;
            }
        }
    }
    // store cell
    if(atomFmtRelative(step.getFmt()) | (write_cell && step.hasCell())){
        auto& cell = j["cell"] = json::object();
        cell["dimension"] = step.getCellDim(step.getFmt() == AtomFmt::Bohr ?
                                            AtomFmt::Bohr :
                                            AtomFmt::Angstrom);
        cell["vectors"] = step.getCellVec();
    }
}

bool JSONWriter(const Molecule &m, std::ostream &file,
                const std::optional<Parameter> &,
                const std::optional<Preset> &p,
                size_t index)
{
    // parse io preset
    if(!p || p->getFmt() != &Plugins::JSON){
        throw IOError("JSON-Writer: needs suitable IO preset");
    }
    auto pp = *p;
    auto indent = std::get<bool>(pp.at("indent").first);
    auto bonds = std::get<bool>(pp.at("bonds").first);
    auto cell = std::get<bool>(pp.at("cell").first);
    auto trajec = std::get<bool>(pp.at("trajectory").first);
    auto fmt = std::get<NamedEnum>(pp.at("fmt").first);
    auto elements = std::get<NamedEnum>(pp.at("elements").first);
    // setup json
    auto j = json{};
    // write periodic table if requested
    if(elements.value() != 0){
        auto pte = m.getPTE();
        if(elements.value() == 1){ // custom
            // remove types not in global PTE
            const auto& root = pte.root ? *pte.root : Vipster::pte;
            for(auto it = pte.begin(); it != pte.end();){
                if(root.find(it->first) != root.end()){
                    it = pte.erase(it);
                }else{
                    ++it;
                }
            }
        }
        else if(elements.value() == 2){ // used
            // remove unused types
            const auto &types = m.getTypes();
            for(auto it = pte.begin(); it != pte.end();){
                if(types.find(it->first) == types.end()){
                    it = pte.erase(it);
                }else{
                    ++it;
                }
            }
        }
        auto &tgt = j["elements"] = json::object();
        for(const auto &[name, elem]: pte){
            tgt[name] = elem;
        }
    }
    // write molecule
    if(trajec){
        auto &mol = j["molecule"] = json::array();
        // iterate over relevant steps
        for (const auto &step: m.getSteps()){
            mol.push_back(json{});
            writeStep(mol.back(), step, fmt, cell, bonds);
        }
    }else{
        auto &mol = j["molecule"] = json{};
        writeStep(mol, m.getStep(index), fmt, cell, bonds);
    }
    file << j.dump(indent ? 0 : -1);
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
