#include "lmpinput.h"
#include "../util.h"

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <sstream>

using namespace Vipster;

// TODO: support velocities?
// TODO: parse coeffs?

static IO::Preset makePreset()
{
    return {&IO::LmpInput,
        {{"style", {NamedEnum{0, {"angle", "atomic", "bond",
                                  "charge", "full", "molecular"}},
                    "The atom_style to be used by Lammps.\n"
                    "See https://lammps.sandia.gov/doc/atom_style.html for details."}},
         {"bonds", {false, "Toggle printing bonds (taken directly from Step)"}},
         {"angles", {false, "Toggle printing angles (extrapolated from bond network)"}},
         {"dihedrals", {false, "Toggle printing dihedrals (extrapolated from bond network)"}},
         {"impropers", {false, "Toggle printing impropers (extrapolated from bond network)"}},
         {"coeff", {false, "Toggle printing parameters from parameter set"}}
    }};
}

using coeffmap = std::map<std::string, std::string>;

static IO::Parameter makeParameter()
{
    return {&IO::LmpInput, {
            {"Pair Coeff", {coeffmap{},
                "If the IOPreset enables parameter printing, "
                "pair coefficients are looked up according the atom type. "
                "Missing types will cause an error."}},
            {"Bond Coeff", {coeffmap{},
                "If the IOPreset enables both bond and parameter printing, "
                "parameters will be looked up according to the bond type. "
                "Missing types will use the fallback (if present) or cause an error."}},
            {"Bond Coeff Fallback", {std::string{},
                "If this value is present, it will be used if a bond type can not "
                "be found in the Bond Coeff map. If it is absent, missing bond types are an error."}},
            {"Angle Coeff", {coeffmap{},
                "If the IOPreset enables both angle and parameter printing, "
                "parameters will be looked up according to the angle type. "
                "Missing types will use the fallback (if present) or cause an error."}},
            {"Angle Coeff Fallback", {std::string{},
                "If this value is present, it will be used if a angle type can not "
                "be found in the Angle Coeff map. If it is absent, missing angle types are an error."}},
            {"Dihedral Coeff", {coeffmap{},
                "If the IOPreset enables both dihedral and parameter printing, "
                "parameters will be looked up according to the dihedral type. "
                "Missing types will use the fallback (if present) or cause an error."}},
            {"Dihedral Coeff Fallback", {std::string{},
                "If this value is present, it will be used if a dihedral type can not "
                "be found in the Dihedral Coeff map. If it is absent, missing dihedral types are an error."}},
            {"Improper Coeff", {coeffmap{},
                "If the IOPreset enables both improper and parameter printing, "
                "parameters will be looked up according to the improper type. "
                "Missing types will use the fallback (if present) or cause an error."}},
            {"Improper Coeff Fallback", {std::string{},
                "If this value is present, it will be used if a improper type can not "
                "be found in the Improper Coeff map. If it is absent, missing improper types are an error."}},
        }};
}

enum class lmpTok{
    type,
    pos,
    charge,
    mol,
    ignore
};

const static std::map<std::string, std::vector<lmpTok>> fmtmap{
    {"angle", {{lmpTok::mol, lmpTok::type, lmpTok::pos}}},
    {"atomic", {{lmpTok::type, lmpTok::pos}}},
    {"body", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"bond", {{lmpTok::mol, lmpTok::type, lmpTok::pos}}},
    {"charge", {{lmpTok::type, lmpTok::charge, lmpTok::pos}}},
    {"dipole", {{lmpTok::type, lmpTok::charge, lmpTok::pos}}},
    {"dpd", {{lmpTok::type, lmpTok::ignore, lmpTok::pos}}},
    {"edpd", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"mdpd", {{lmpTok::type, lmpTok::ignore, lmpTok::pos}}},
    {"tdpd", {{lmpTok::type, lmpTok::pos}}},
    {"electron", {{lmpTok::type, lmpTok::charge, lmpTok::ignore,
                   lmpTok::ignore, lmpTok::pos}}},
    {"ellipsoid", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"full", {{lmpTok::mol, lmpTok::type, lmpTok::charge, lmpTok::pos}}},
    {"line", {{lmpTok::mol, lmpTok::type, lmpTok::ignore,
               lmpTok::ignore, lmpTok::pos}}},
    {"meso", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore,
               lmpTok::ignore, lmpTok::pos}}},
    {"molecular", {{lmpTok::mol, lmpTok::type, lmpTok::pos}}},
    {"peri", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"smd", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::ignore,
              lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"sphere", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"template", {{lmpTok::mol, lmpTok::ignore, lmpTok::ignore,
                   lmpTok::type, lmpTok::pos}}},
    {"tri", {{lmpTok::mol, lmpTok::type, lmpTok::ignore,
              lmpTok::ignore, lmpTok::pos}}},
    {"wavepacket", {{lmpTok::type, lmpTok::charge, lmpTok::ignore,lmpTok::ignore,
                     lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"hybrid", {{lmpTok::type, lmpTok::pos}}}
};

std::vector<lmpTok> getFmtGuess(std::istream& file, size_t nat){
    // BUG: WILL fail if fmt == tdpd, hybrid, template
    // probably also for dipole and ellipsoid
    auto rewindpos = file.tellg();
    std::vector<std::vector<std::string>> atoms;
    std::string line, tok;
    // store tokenized atom-lines
    nat = std::min(nat, static_cast<size_t>(20));
    atoms.resize(nat);
    for(auto& at: atoms) {
        std::getline(file, line);
        std::stringstream ss{line};
        while((ss >> tok)){
            if(tok[0] == '#') break;
            at.push_back(tok);
        }
    }
    // determine base argument length
    size_t narg{100};
    for(const auto& at:atoms){
        narg = std::min(narg,at.size());
    }
    file.seekg(rewindpos);
    // collect parser-pieces
    auto checkInt = [&atoms](size_t col){
        try{
            for (auto& at: atoms){
                if(at[col].find('.') != at[col].npos){
                    return false;
                }
            }
            return true;
        }catch(...){
            return false;
        }
    };
    if(narg == 5){
        // only one possible fmt (atomic)
        return fmtmap.at("atomic");
    }
    if(narg == 6){
        if(checkInt(2)){
            // angle/molecular have molID, then type
            return fmtmap.at("angle");
        }
        // parse as charge, even though this is most likely wrong
        return fmtmap.at("charge");
    }
    std::vector<lmpTok> parser{};
    size_t col{1}, poscoord{narg-4};
    /* assume:
     * - trailing int-columns are image-flags
     * - three cols before image-flags are position
     * - second or third col are atomtype
     * - first col between type and pos is charge (if present)
     */
    if (checkInt(2)) {
        // assume col1 is molID, probably fails for ellipsoid
        parser.push_back(lmpTok::ignore);
        col++;
    }
    parser.push_back(lmpTok::type);
    for(size_t img=narg-1; img >= std::max(narg-3, static_cast<size_t>(5)) &&
                           poscoord-col >= 0 ; --img){
        if (!checkInt(img)){
            break;
        }
        --poscoord;
    }
    if((poscoord-col)>0){
        parser.push_back(lmpTok::charge);
        ++col;
        for(size_t i=0; i<(poscoord-col); ++i){
            parser.push_back(lmpTok::ignore);
        }
    }
    parser.push_back(lmpTok::pos);
    return parser;
}

auto makeParser(std::vector<lmpTok> fmt){
    return [fmt](std::istream& file, Step& s, size_t nat,
                 std::map<size_t, std::string>& types,
                 std::map<size_t, size_t>& indices){
        s.newAtoms(nat);
        std::string line{}, dummy{};
        size_t id{};
        for (auto it=s.begin(); it!=s.end(); ++it) {
            auto& at = *it;
            std::getline(file, line);
            std::stringstream ss{line};
            ss >> id;
            indices[id] = it.getIdx();
            for (lmpTok tok: fmt) {
                switch(tok){
                case lmpTok::type:
                    ss >> id;
                    at.name = types[id];
                    break;
                case lmpTok::charge:
                    ss >> at.properties->charge;
                    break;
                case lmpTok::pos:
                    ss >> at.coord[0] >> at.coord[1] >> at.coord[2];
                    break;
                case lmpTok::mol:
                    [[fallthrough]];
                case lmpTok::ignore:
                    ss >> dummy;
                    break;
                }
            }
        }
    };
}

void groupSets(std::list<std::set<size_t>>& molecules){
    auto size = molecules.size();
    std::set<size_t> test;
    // pairwise compare sets
    for(auto it1 = molecules.begin(); it1 != molecules.end(); ++it1){
        auto it2 = it1;
        ++it2;
        while(it2 != molecules.end()){
            test.clear();
            std::set_intersection(it1->begin(), it1->end(),
                                  it2->begin(), it2->end(),
                                  std::inserter(test, test.begin()));
            if(!test.empty()){
                it1->insert(it2->begin(), it2->end());
                it2 = molecules.erase(it2);
            }else{
                ++it2;
            }
        }
    }
    if(molecules.size() != size){
        groupSets(molecules);
    }
};

auto makeWriter(const std::vector<lmpTok>& fmt,
                const std::vector<size_t>& molID,
                const std::map<std::string, size_t>& atomtypemap){
    return [&](std::ostream& file, const auto& step){
        for(auto it=step.begin(); it!=step.end(); ++it){
            file << std::left << std::setw(3) << (it.getIdx()+1);
            for(const auto& tok: fmt){
                file << ' ';
                switch(tok){
                case lmpTok::charge:
                    file << std::left << std::setw(3)
                         << (*it).properties->charge;
                    break;
                case lmpTok::mol:
                    file << std::left << std::setw(3)
                         << molID[it.getIdx()];
                    break;
                case lmpTok::type:
                    file << std::left << std::setw(3)
                         << atomtypemap.at((*it).name);
                    break;
                case lmpTok::pos:
                    file << std::right
                         << std::setw(std::numeric_limits<Vec::value_type>::max_digits10+5) << (*it).coord[0]
                         << ' ' << std::setw(std::numeric_limits<Vec::value_type>::max_digits10+5) << (*it).coord[1]
                         << ' ' << std::setw(std::numeric_limits<Vec::value_type>::max_digits10+5) << (*it).coord[2];
                    break;
                case lmpTok::ignore:
                    break;
                }
            }
            file << '\n' << std::left;
        }
    };
}

IO::Data LmpInpParser(const std::string& name, std::istream &file)
{
    enum class ParseMode{Header,Atoms,Types};

    IO::Data data{};
    Molecule& m = data.mol;
    m.name = name;
    Step& s = m.newStep();
    s.setFmt(AtomFmt::Angstrom);
    s.setCellDim(1, AtomFmt::Angstrom);

    std::string tmp;
    size_t nat{}, nbnd{}, ntype{};
    double t1, t2;
    Mat cell{};
    std::map<size_t, size_t> indices{};
    std::map<size_t, std::string> types{};
    std::map<size_t, std::string> bondtypes{};
    while (std::getline(file, tmp)) {
        auto [line, comment] = IO::stripComment(tmp);
        if (line.find("atoms") != std::string::npos) {
            std::stringstream ss{line};
            ss >> nat;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse number of atoms");
            }
        } else if (line.find("atom types") != std::string::npos) {
            std::stringstream ss{line};
            ss >> ntype;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse number of types");
            }
        } else if (line.find("bonds") != std::string::npos) {
            // bonds given -> deactive automatic bond generation
            std::stringstream ss{line};
            ss >> nbnd;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse number of bonds");
            }
        } else if (line.find("bond types") != std::string::npos){
            // try comment-only lines to find out named bond types
            bool seekType{true};
            auto rewindpos = file.tellg();
            while(seekType && std::getline(file, tmp)){
                std::tie(line, comment) = IO::stripComment(tmp);
                if(line.empty() && !comment.empty()){
                    // use first value as idx, rest of line as name
                    size_t pos{};
                    size_t idx = std::stoul(comment, &pos);
                    bondtypes[idx] = IO::trim(comment.substr(pos));
                    // if succesfull, we consume this line
                    rewindpos = file.tellg();
                }else{
                    seekType = false;
                }
            }
            file.seekg(rewindpos);
        } else if (line.find("xlo xhi") != std::string::npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse cell X dimension");
            }
            cell[0][0] = t2 - t1;
        } else if (line.find("ylo yhi") != std::string::npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse cell Y dimension");
            }
            cell[1][1] = t2 - t1;
        } else if (line.find("zlo zhi") != std::string::npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse cell Z dimension");
            }
            cell[2][2] = t2 - t1;
        } else if (line.find("xy xz yz") != std::string::npos) {
            std::stringstream ss{line};
            ss >> cell[1][0] >> cell[2][0] >> cell[2][1];
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to "
                                "parse cell tilt factors");
            }
        } else if (line.find("Masses") != std::string::npos) {
            // skip empty line
            std::getline(file, line);
            size_t id;
            std::string name;
            for (size_t i=0; i<ntype; ++i) {
                std::getline(file, line);
                std::stringstream ss{line};
                ss >> id >> t1;
                std::size_t cpos = line.find('#');
                if(cpos != std::string::npos) {
                    // if there's a comment, extract the typename from it
                    std::stringstream ss{line.substr(cpos+1)};
                    ss >> name;
                    types[id] = name;
                } else {
                    // else just number the types accordingly
                    name = std::to_string(id);
                    types[id] = name;
                    // and try to guess type from mass
                    m.getPTE().insert_or_assign(name,
                        [&t1](){
                        const Vipster::PeriodicTable::mapped_type* cur_guess{&Vipster::pte.at("")};
                        double best_diff{5};
                        for(const auto& [name, type]: Vipster::pte){
                            double cur_diff = std::abs(t1-type.m);
                            if(cur_diff < best_diff){
                                best_diff = cur_diff;
                                cur_guess = &type;
                            }
                        }
                        return *cur_guess;
                    }());
                }
                if (ss.fail()) {
                    throw IO::Error("Lammps Input: failed to parse atom type");
                }
                // finally, save mass
                m.getPTE()[types[id]].m = t1;
            }
        } else if (line.find("Atoms") != std::string::npos) {
            std::vector<lmpTok> fmt{};
            // lookup fixed parser if format is given
            if (!comment.empty()) {
                auto pos = fmtmap.find(IO::trim(comment));
                if(pos != fmtmap.end()){
                    fmt = pos->second;
                }
            }
            // skip empty line
            std::getline(file, line);
            // if no format was given, try to determine a suitable parser
            if (fmt.empty()) {
                fmt = getFmtGuess(file, nat);
            }
            // do the parsing
            makeParser(fmt)(file, s, nat, types, indices);
            // cell needs to be known at this point, so let's set it
            s.setCellVec(cell);
        } else if (line.find("Bonds") != std::string::npos) {
            // skip empty line
            std::getline(file, line);
            // operate on crystal coordinates to interpret periodic bonds more easily
            auto sc = s.asFmt(AtomFmt::Crystal);
            for(size_t i=0; i<nbnd; ++i){
                size_t n, id, at1, at2;
                std::getline(file, line);
                std::stringstream{line} >> n >> id >> at1 >> at2;
                // lookup atoms by file-local index
                auto f1 = indices.find(at1);
                auto f2 = indices.find(at2);
                if ((f1 == indices.end()) || (f2 == indices.end())) {
                    throw IO::Error{"Lammps Input: invalid atom IDs in bond: "+line};
                }
                at1 = f1->second;
                at2 = f2->second;
                // check if bond should be periodic, save accordingly
                auto dist = sc.at(at1).coord - sc.at(at2).coord;
                DiffVec diff, dir;
                // diff contains integer distance in cell-units
                std::transform(dist.begin(), dist.end(), diff.begin(), truncf);
                // dist contains distance inside of cell
                std::transform(dist.begin(), dist.end(), dist.begin(),
                    [](double f){return std::fmod(f,1);});
                // dir contains direction of dist
                std::transform(dist.begin(), dist.end(), dir.begin(),
                    [](double f){
                        return (std::abs(f) < std::numeric_limits<double>::epsilon())?
                                    0 : ((f<0) ? -1 : 1);
                    });
                // fail if atoms overlap
                if(std::none_of(dir.begin(), dir.end(), [](auto& i){return i!=0;}) &&
                   std::none_of(diff.begin(), diff.end(), [](auto& i){return i!=0;})){
                    char errmsg[50];
                    sprintf(errmsg, "Lammps Input: failed to parse bond %lu", i);
                    throw IO::Error{errmsg};
                }
                // wrap if needed
                for(size_t d=0; d<3; ++d){
                    if(std::abs(dist[d]) > 0.5){
                        dist[d] -= dir[d];
                        diff[d] += dir[d];
                    }
                }
                // if bond type has been annotated, use name
                auto pos = bondtypes.find(id);
                if(pos != bondtypes.end()){
                    sc.addBond(at1, at2, diff, pos->second);
                }else{
                    sc.addBond(at1, at2, diff, std::to_string(id));
                }
            }
        }
    }
    return data;
}

bool LmpInpWriter(const Molecule& m, std::ostream &file,
                  const std::optional<IO::Parameter>& p,
                  const std::optional<IO::Preset>& c,
                  size_t index)
{
    const auto step = m.getStep(index).asFmt(AtomFmt::Angstrom);
    // parse iopreset
    if(!c || c->getFmt() != &IO::LmpInput){
        throw IO::Error("Lammps Input: writer needs suitable IO preset");
    }
    auto print_coeff = std::get<bool>(c->at("coeff").first);
    auto print_bonds = std::get<bool>(c->at("bonds").first);
    auto print_angles = std::get<bool>(c->at("angles").first);
    auto print_dihed = std::get<bool>(c->at("dihedrals").first);
    auto print_improp = std::get<bool>(c->at("impropers").first);
    auto [angles, dihedrals, impropers] = step.getTopology(print_angles, print_dihed, print_improp);
    const auto &style = std::get<NamedEnum>(c->at("style").first);
    const auto tokens = fmtmap.at(style.at(style));
    bool needsMolID = std::find(tokens.begin(), tokens.end(), lmpTok::mol) != tokens.end();

    file << std::setprecision(std::numeric_limits<Vec::value_type>::max_digits10);

    /* prepare topology
     *
     * *typemap contains names and their lammps-internal type index, starting with 1
     * *typelist maps the lammps type index to the vipster item index
     */
    std::vector<size_t> bondtypelist;
    std::map<std::string, size_t> bondtypemap;
    if(print_bonds){
        for(const auto& bond: step.getBonds()){
            if(bond.type){
                bondtypelist.push_back(bondtypemap.emplace(
                    bond.type->first,
                    bondtypemap.size()+1).first->second);
            }else{
                const std::string& name1 = step[bond.at1].name;
                const std::string& name2 = step[bond.at2].name;
                bondtypelist.push_back(bondtypemap.emplace(
                    fmt::format("{}-{}", std::min(name1, name2), std::max(name1, name2)),
                    bondtypemap.size()+1).first->second);
            }
        }
    }

    std::map<std::string, size_t> angletypemap;
    std::vector<size_t> angletypelist;
    if(print_angles){
        for(const auto& angle: angles){
            const std::string& name1 = step[angle.at1].name;
            const std::string& name2 = step[angle.at2].name;
            const std::string& name3 = step[angle.at3].name;
            angletypelist.push_back(angletypemap.emplace(
                fmt::format("{}-{}-{}", std::min(name1, name3), name2, std::max(name1, name3)),
                angletypemap.size()+1).first->second);
        }
    }

    std::vector<size_t> improptypelist;
    std::map<std::string, size_t> improptypemap;
    if(print_improp){
        for(const auto& improp: impropers){
            const std::string& name1 = step[improp.at1].name.c_str();
            std::array<std::string_view, 3> names = {step[improp.at2].name,
                                                     step[improp.at3].name,
                                                     step[improp.at4].name};
            std::sort(names.begin(), names.end());
            improptypelist.push_back(improptypemap.emplace(
                fmt::format("{}-{}-{}-{}", name1, names[0], names[1], names[2]),
                improptypemap.size()+1).first->second);
        }
    }

    std::vector<size_t> dihedtypelist;
    std::map<std::string, size_t> dihedtypemap; // map vipster type to lammps type
    if(print_dihed){
        for(const auto& dihed: dihedrals){
            const std::string& name1 = step[dihed.at1].name;
            const std::string& name2 = step[dihed.at2].name;
            const std::string& name3 = step[dihed.at3].name;
            const std::string& name4 = step[dihed.at4].name;
            dihedtypelist.push_back(dihedtypemap.emplace(
                name1 < name4 ?
                    fmt::format("{}-{}-{}-{}", name1, name2, name3, name4) :
                    fmt::format("{}-{}-{}-{}", name4, name3, name2, name1),
                dihedtypemap.size()+1).first->second);
        }
    }

    // prepare Molecule-IDs
    std::vector<size_t> molID(step.getNat());
    std::list<std::set<size_t>> molSets{};
    if(needsMolID){
        molID.resize(step.getNat());
        for(const auto& bond: step.getBonds()){
            molSets.push_back(std::set<size_t>{bond.at1, bond.at2});
        }
        groupSets(molSets);
        auto it = molSets.begin();
        for(size_t i=0; i<molSets.size(); ++i, ++it){
            for(const auto& at: *it){
                molID[at] = i+1;
            }
        }
    }

    /*
     * Header
     */
    file << '\n'
         << step.getNat() << " atoms\n"
         << step.getNtyp() << " atom types\n";
    if(print_bonds && !bondtypelist.empty()){
        file << bondtypelist.size() << " bonds\n"
             << bondtypemap.size() << " bond types\n";
        for(const auto& pair: bondtypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }
    if(print_angles && !angletypelist.empty()){
        file << angletypelist.size() << " angles\n"
             << angletypemap.size() << " angle types\n";
        for(const auto& pair: angletypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }
    if(print_dihed && !dihedtypelist.empty()){
        file << dihedtypelist.size() << " dihedrals\n"
             << dihedtypemap.size() << " dihedral types\n";
        for(const auto& pair: dihedtypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }
    if(print_improp && !improptypelist.empty()){
        file << improptypelist.size() << " impropers\n"
             << improptypemap.size() << " improper types\n";
        for(const auto& pair: improptypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }

    if(step.hasCell()){
        // print existing cell
        auto vec = step.getCellVec() * step.getCellDim(AtomFmt::Angstrom);
        if(!float_comp(vec[0][1], 0.) || !float_comp(vec[0][2], 0.) || !float_comp(vec[1][2], 0.)){
            throw IO::Error("Lammps Input: cell vectors must form diagonal or lower triangular matrix");
        }
        file << "\n0.0 "
             << vec[0][0] << " xlo xhi\n0.0 "
             << vec[1][1] << " ylo yhi\n0.0 "
             << vec[2][2] << " zlo zhi\n";
        if(!float_comp(vec[1][0], 0.) || !float_comp(vec[2][0], 0.) || !float_comp(vec[2][1], 0.)){
            file << vec[1][0] << ' ' << vec[2][0] << ' ' << vec[2][1] << " xy xz yz\n";
        }
    }else{
        // create bounding box suitable for shrink-wrapped non-periodic simulations
        Vec pos_min{{std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max()}};
        Vec pos_max{{std::numeric_limits<double>::lowest(),
                 std::numeric_limits<double>::lowest(),
                 std::numeric_limits<double>::lowest()}};
        for(const auto& at: step){
            pos_min[0] = std::min(pos_min[0], at.coord[0]);
            pos_min[1] = std::min(pos_min[1], at.coord[1]);
            pos_min[2] = std::min(pos_min[2], at.coord[2]);
            pos_max[0] = std::max(pos_max[0], at.coord[0]);
            pos_max[1] = std::max(pos_max[1], at.coord[1]);
            pos_max[2] = std::max(pos_max[2], at.coord[2]);
        }
        file << '\n' << std::setprecision(8)
             << pos_min[0] << ' ' << pos_max[0] << " xlo xhi\n"
             << pos_min[1] << ' ' << pos_max[1] << " ylo yhi\n"
             << pos_min[2] << ' ' << pos_max[2] << " zlo zhi\n";
    }

    file << "\nMasses\n\n";
    std::map<std::string, size_t> atomtypemap;
    for(const auto &t: step.getTypes()){
        atomtypemap.emplace(t, atomtypemap.size()+1);
        file << atomtypemap.size() << ' ' << m.getPTE().at(t).m << " # " << t << '\n';
    }

    file << "\nAtoms # " << style.at(style) << "\n\n";
    makeWriter(tokens, molID, atomtypemap)(file, step);

    if(print_bonds && !bondtypelist.empty()){
        file << "\nBonds\n\n";
        const auto& bonds = step.getBonds();
        for(size_t i=0; i!=bondtypelist.size(); ++i){
            file << (i+1) << ' '
                 << bondtypelist[i] << ' '
                 << bonds[i].at1+1 << ' '
                 << bonds[i].at2+1 << '\n';
        }
    }

    if(print_angles && !angletypelist.empty()){
        file << "\nAngles\n\n";
        for(size_t i=0; i!=angletypelist.size(); ++i){
            file << (i+1) << ' '
                 << angletypelist[i] << ' '
                 << angles[i].at1+1 << ' '
                 << angles[i].at2+1 << ' '
                 << angles[i].at3+1 << '\n';
        }
    }

    if(print_dihed && !dihedtypelist.empty()){
        file << "\nDihedrals\n\n";
        for(size_t i=0; i!=dihedtypelist.size(); ++i){
            file << (i+1) << ' '
                 << dihedtypelist[i] << ' '
                 << dihedrals[i].at1+1 << ' '
                 << dihedrals[i].at2+1 << ' '
                 << dihedrals[i].at3+1 << ' '
                 << dihedrals[i].at4+1 << '\n';
        }
    }

    if(print_improp && !improptypelist.empty()){
        file << "\nImpropers\n\n";
        for(size_t i=0; i!=improptypelist.size(); ++i){
            file << (i+1) << ' '
                 << improptypelist[i] << ' '
                 << impropers[i].at1+1 << ' '
                 << impropers[i].at2+1 << ' '
                 << impropers[i].at3+1 << ' '
                 << impropers[i].at4+1 << '\n';
        }
    }

    if(print_coeff){
        // parse parameter
        if(!p || p->getFmt() != &IO::LmpInput){
            throw IO::Error("Lammps Input: print coeffs requested but no suitable parameter set provided");
        }
        const auto& paircoeffs = std::get<coeffmap>(p->at("Pair Coeff").first);
        const auto& bondcoeffs = std::get<coeffmap>(p->at("Bond Coeff").first);
        const auto& bondfallback = IO::trim(std::get<std::string>(p->at("Bond Coeff Fallback").first));
        bool bond_hasfallback = !bondfallback.empty();
        const auto& anglecoeffs = std::get<coeffmap>(p->at("Angle Coeff").first);
        const auto& anglefallback = IO::trim(std::get<std::string>(p->at("Angle Coeff Fallback").first));
        bool angle_hasfallback = !anglefallback.empty();
        const auto& dihedralcoeffs = std::get<coeffmap>(p->at("Dihedral Coeff").first);
        const auto& dihedralfallback = IO::trim(std::get<std::string>(p->at("Dihedral Coeff Fallback").first));
        bool dihedral_hasfallback = !bondfallback.empty();
        const auto& impropercoeffs = std::get<coeffmap>(p->at("Improper Coeff").first);
        const auto& improperfallback = IO::trim(std::get<std::string>(p->at("Improper Coeff Fallback").first));
        bool improper_hasfallback = !bondfallback.empty();

        // Pair Coeffs
        file << "\nPair Coeffs\n\n";
        // use atomtypemap created in the Masses section
        for(const auto &[t, idx]: atomtypemap){
            auto pos = paircoeffs.find(t);
            if(pos == paircoeffs.end()){
                throw IO::Error{fmt::format("Could not find pair coefficients for atom type {} in parameter set.", t)};
            }
            fmt::print(file, "{} {}\n", idx, pos->second);
        }

        // Bond Coeffs
        if(print_bonds){
            file << "\nBond Coeffs\n\n";
            for(const auto &[t, idx]: bondtypemap){
                auto pos = bondcoeffs.find(t);
                if(pos == bondcoeffs.end()){
                    if(bond_hasfallback){
                        fmt::print(file, "{} {}\n", idx, bondfallback);
                    }else{
                        throw IO::Error{fmt::format("Could not find bond coefficients for bond type {} in parameter set. "
                                                    "No fallback value provided", t)};
                    }
                }else{
                    fmt::print(file, "{} {}\n", idx, pos->second);
                }
            }
        }

        // Angle Coeffs
        if(print_angles){
            file << "\nAngle Coeffs\n\n";
            for(const auto &[t, idx]: angletypemap){
                auto pos = anglecoeffs.find(t);
                if(pos == anglecoeffs.end()){
                    if(angle_hasfallback){
                        fmt::print(file, "{} {}\n", idx, anglefallback);
                    }else{
                        throw IO::Error{fmt::format("Could not find angle coefficients for angle type {} in parameter set. "
                                                    "No fallback value provided", t)};
                    }
                }else{
                    fmt::print(file, "{} {}\n", idx, pos->second);
                }
            }
        }

        // Dihedral Coeffs
        if(print_dihed){
            file << "\nDihedral Coeffs\n\n";
            for(const auto &[t, idx]: dihedtypemap){
                auto pos = dihedralcoeffs.find(t);
                if(pos == dihedralcoeffs.end()){
                    if(dihedral_hasfallback){
                        fmt::print(file, "{} {}\n", idx, dihedralfallback);
                    }else{
                        throw IO::Error{fmt::format("Could not find dihedral coefficients for dihedral type {} in parameter set. "
                                                    "No fallback value provided", t)};
                    }
                }else{
                    fmt::print(file, "{} {}\n", idx, pos->second);
                }
            }
        }

        // Improper Coeffs
        if(print_improp){
            file << "\nImproper Coeffs\n\n";
            for(const auto &[t, idx]: improptypemap){
                auto pos = impropercoeffs.find(t);
                if(pos == impropercoeffs.end()){
                    if(improper_hasfallback){
                        fmt::print(file, "{} {}\n", idx, improperfallback);
                    }else{
                        throw IO::Error{fmt::format("Could not find improper coefficients for improper type {} in parameter set. "
                                                    "No fallback value provided", t)};
                    }
                }else{
                    fmt::print(file, "{} {}\n", idx, pos->second);
                }
            }
        }
    }

    file << '\n';
    return true;
}

const IO::Plugin IO::LmpInput =
{
    "Lammps Data File",
    "lmp",
    "lmp",
    &LmpInpParser,
    &LmpInpWriter,
    &makeParameter,
    &makePreset
};
