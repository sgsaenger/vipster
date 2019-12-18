#include "plugin.h"
#include "../util.h"
#include <sstream>

using namespace Vipster;

enum class LMPAtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};

enum class lmpTok{
    type,
    pos,
    charge,
    mol,
    ignore
};

const static std::map<LMPAtomStyle, std::string> fmt2str{
    {LMPAtomStyle::Angle, "angle"},
    {LMPAtomStyle::Atomic, "atomic"},
    {LMPAtomStyle::Bond, "bond"},
    {LMPAtomStyle::Charge, "charge"},
    {LMPAtomStyle::Full, "full"},
    {LMPAtomStyle::Molecular, "molecular"},
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
    // TODO: WILL fail if fmt == tdpd, hybrid, template
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
                auto f = stod(at[col]);
                if(f != floor(f)){
                    return false;
                }
            }
            return true;
        }catch(...){
            return false;
        }
    };
    if(narg == 5){
        // only one possible setup (atomic)
        return {lmpTok::type, lmpTok::pos};
    }
    if(narg == 6){
        if(checkInt(2)){
            // angle/molecular have molID, then type
            return {lmpTok::ignore, lmpTok::type, lmpTok::pos};
        }
        // parse as charge, even though this is most likely wrong
        return {lmpTok::type, lmpTok::charge, lmpTok::pos};
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
    for(size_t img=narg-1; img>=std::max(narg-3,static_cast<size_t>(5)); --img){
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
    return [&](std::ostream& file, const Step& s){
        for(auto it=s.begin(); it!=s.end(); ++it){
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
    m.setName(name);
    Step& s = m.newStep();
    s.setFmt(AtomFmt::Angstrom);
    s.setCellDim(1, CdmFmt::Angstrom);

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
            s.setBondMode(BondMode::Manual);
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
                    s.pte->insert_or_assign(name,
                        [&t1](){
                        const Vipster::PeriodicTable::mapped_type* cur_guess{&Vipster::pte.at("")};
                        double best_diff{5};
                        for(const auto& pair: Vipster::pte){
                            double cur_diff = std::abs(t1-pair.second.m);
                            if(cur_diff < best_diff){
                                best_diff = cur_diff;
                                cur_guess = &pair.second;
                            }
                        }
                        return *cur_guess;
                    }());
                }
                if (ss.fail()) {
                    throw IO::Error("Lammps Input: failed to parse atom type");
                }
                // finally, save mass
                (*s.pte)[types[id]].m = t1;
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
                    throw IO::Error{"Invalid atom IDs in bond: "+line};
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
                    sc.newBond(at1, at2, diff, pos->second);
                }else{
                    sc.newBond(at1, at2, diff, std::to_string(id));
                }
            }
        }
    }
    return data;
}

bool LmpInpWriter(const Molecule& m, std::ostream &file,
                  const IO::BaseParam *const,
                  const std::optional<IO::BasePreset>& c,
                  size_t index)
{
    const auto step = m.getStep(index).asFmt(AtomFmt::Angstrom);
    // parse iopreset
    if(!c || c->getFmt() != &IO::LmpInput){
        throw IO::Error("Lammps-Writer needs suitable IO preset");
    }
    auto bonds = std::get<bool>(c->at("bonds"));
    auto angles = std::get<bool>(c->at("angles"));
    auto dihedrals = std::get<bool>(c->at("dihedrals"));
    auto impropers = std::get<bool>(c->at("impropers"));
    auto style = static_cast<LMPAtomStyle>(std::get<uint>(c->at("style")));
    const auto tokens = fmtmap.at(fmt2str.at(style));
    bool needsMolID = std::find(tokens.begin(), tokens.end(), lmpTok::mol) != tokens.end();

    file << std::setprecision(std::numeric_limits<Vec::value_type>::max_digits10);

    // prepare bonds
    std::vector<std::tuple<size_t, size_t, size_t>> bondlist;
    std::map<std::string, size_t> bondtypemap;
    if(bonds || angles || dihedrals || impropers || needsMolID){
        for(auto& bond: step.getBonds()){
            if(bond.type){
                bondlist.push_back({std::min(bond.at1, bond.at2),
                                    std::max(bond.at1, bond.at2),
                                    bondtypemap.emplace(bond.type->first,
                                     bondtypemap.size()+1).first->second});
            }else{
                const std::string& t1 = step[bond.at1].name;
                const std::string& t2 = step[bond.at2].name;
                bondlist.push_back({std::min(bond.at1, bond.at2),
                                    std::max(bond.at1, bond.at2),
                                    bondtypemap.emplace(
                                     std::min(t1, t2)+'-'+std::max(t1,t2),
                                     bondtypemap.size()+1).first->second});
            }
        }
        std::sort(bondlist.begin(), bondlist.end(), [](const auto& l, const auto& r){
            return std::tie(std::get<0>(l), std::get<1>(l)) <
                   std::tie(std::get<0>(r), std::get<1>(r));
        });
    }

    // prepare angles and impropers
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> anglelist;
    std::map<std::string, size_t> angletypemap;
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> improplist;
    std::vector<std::string> improptypelist;
    std::map<std::string, size_t> improptypemap;
    if(angles || dihedrals || impropers){
        for(auto it=bondlist.begin(); it!=bondlist.end(); ++it){
            const auto& a0 = std::get<0>(*it);
            const auto& a1 = std::get<1>(*it);
            const std::string& na0 = step[a0].name;
            const std::string& na1 = step[a1].name;
            for(auto it2=it+1; it2!=bondlist.end(); ++it2){
                const auto& b0 = std::get<0>(*it2);
                const auto& b1 = std::get<1>(*it2);
                const std::string& nb0 = step[b0].name;
                const std::string& nb1 = step[b1].name;
                bool found{false};
                if(a0 == b0){
                    anglelist.emplace_back(a1, a0, b1,
                        angletypemap.emplace(
                          std::min(na1, nb1)+'-'+na0+'-'+std::max(na1, nb1),
                          angletypemap.size()+1).first->second);
                    found = true;
                }else if(a1 == b0){
                    anglelist.emplace_back(a0, a1, b1,
                        angletypemap.emplace(
                          std::min(na0, nb1)+'-'+na1+'-'+std::max(na0, nb1),
                          angletypemap.size()+1).first->second);
                    found = true;
                }else if(a1 == b1){
                    anglelist.emplace_back(a0, a1, b0,
                        angletypemap.emplace(
                          std::min(na0, nb0)+'-'+na1+'-'+std::max(na0, nb0),
                          angletypemap.size()+1).first->second);
                    found = true;
                }
                // a0 == b1 impossible because of sorting
                if(impropers && found){
                    for(auto it3=it2+1; it3!=bondlist.end(); ++it3){
                        const auto& c0 = std::get<0>(*it3);
                        const auto& c1 = std::get<1>(*it3);
                        const auto& ang = anglelist.back();
                        const auto& t0 = std::get<0>(ang);
                        const auto& t1 = std::get<1>(ang);
                        const auto& t2 = std::get<2>(ang);
                        if(c0 == t1){
                            improplist.emplace_back(t1, t0, t2, c1);
                        }else if(c1 == t1){
                            improplist.emplace_back(t1, t0, t2, c0);
                        }
                    }
                }
            }
        }
        for(const auto& improp: improplist){
            std::string t0 = step[std::get<0>(improp)].name;
            std::string t1 = step[std::get<1>(improp)].name;
            std::string t2 = step[std::get<2>(improp)].name;
            std::string t3 = step[std::get<3>(improp)].name;
            if(t1 < t2){
                std::swap(t1, t2);
            }
            if(t1 < t3){
                std::swap(t1, t3);
            }
            if(t2 < t3){
                std::swap(t2, t3);
            }
            improptypelist.emplace_back(t0+'-'+t1+'-'+t2+'-'+t3);
            improptypemap.emplace(improptypelist.back(), improptypemap.size()+1);
        }
    }

    // prepare dihedrals
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> dihedlist;
    std::vector<std::string> dihedtypelist;
    std::map<std::string, size_t> dihedtypemap;
    if(dihedrals){
        for(auto it=anglelist.begin(); it!=anglelist.end(); ++it){
            const auto& a0 = std::get<0>(*it);
            const auto& a1 = std::get<1>(*it);
            const auto& a2 = std::get<2>(*it);
            for(auto it2=it+1; it2!=anglelist.end(); ++it2){
                const auto& b0 = std::get<0>(*it2);
                const auto& b1 = std::get<1>(*it2);
                const auto& b2 = std::get<2>(*it2);
                if(a1 == b0){
                    if(b1 == a0){
                        if(a2 != b2){
                            dihedlist.emplace_back(a2, a1, a0, b2);
                        }
                    }else if(b1 == a2){
                        if(a0 != b2){
                            dihedlist.emplace_back(a0, a1, a2, b2);
                        }
                    }
                }else if(a1 == b2){
                    if(b1 == a0){
                        if(a2 != b0){
                            dihedlist.emplace_back(a2, a1, a0, b0);
                        }
                    }else if(b1 == a2){
                        if(a0 != b0){
                            dihedlist.emplace_back(a0, a1, a2, b0);
                        }
                    }
                }
            }
        }
        for(const auto& dihed: dihedlist){
            const std::string& t0 = step[std::get<0>(dihed)].name;
            const std::string& t1 = step[std::get<1>(dihed)].name;
            const std::string& t2 = step[std::get<2>(dihed)].name;
            const std::string& t3 = step[std::get<3>(dihed)].name;
            if((t0 < t3) || ((t0 == t3) && (t1 < t2))){
                dihedtypelist.emplace_back(t0+'-'+t1+'-'+t2+'-'+t3);
            }else{
                dihedtypelist.emplace_back(t3+'-'+t2+'-'+t1+'-'+t0);
            }
            dihedtypemap.emplace(dihedtypelist.back(), dihedtypemap.size()+1);
        }
    }

    // prepare Molecule-IDs
    std::vector<size_t> molID(step.getNat());
    std::list<std::set<size_t>> molSets{};
    if(needsMolID){
        molID.resize(step.getNat());
        // make sure this doesn't reset bonds so the user won't be surprised
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
    if(bonds && !bondlist.empty()){
        file << bondlist.size() << " bonds\n"
             << bondtypemap.size() << " bond types\n";
        for(const auto& pair: bondtypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }
    if(angles && !anglelist.empty()){
        file << anglelist.size() << " angles\n"
             << angletypemap.size() << " angle types\n";
        for(const auto& pair: angletypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }
    if(dihedrals && !dihedlist.empty()){
        file << dihedlist.size() << " dihedrals\n"
             << dihedtypemap.size() << " dihedral types\n";
        for(const auto& pair: dihedtypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }
    if(impropers && !improplist.empty()){
        file << improplist.size() << " impropers\n"
             << improptypemap.size() << " improper types\n";
        for(const auto& pair: improptypemap){
            file << '#' << pair.second << ' ' << pair.first << '\n';
        }
    }

    auto vec = step.getCellVec() * step.getCellDim(CdmFmt::Angstrom);
    if(!float_comp(vec[0][1], 0.) || !float_comp(vec[0][2], 0.) || !float_comp(vec[1][2], 0.)){
        throw IO::Error("Cell vectors must form diagonal or lower triangular matrix for Lammps");
    }
    file << "\n0.0 "
         << vec[0][0] << " xlo xhi\n0.0 "
         << vec[1][1] << " ylo yhi\n0.0 "
         << vec[2][2] << " zlo zhi\n";
    if(!float_comp(vec[1][0], 0.) || !float_comp(vec[2][0], 0.) || !float_comp(vec[2][1], 0.)){
        file << vec[1][0] << ' ' << vec[2][0] << ' ' << vec[2][1] << " xy xz yz\n";
    }

    file << "\nMasses\n\n";
    std::map<std::string, size_t> atomtypemap;
    for(const auto& t: step.getTypes()){
        atomtypemap.emplace(t, atomtypemap.size()+1);
        file << atomtypemap.size() << ' ' << step.pte->at(t).m << " # " << t << '\n';
    }

    file << "\nAtoms # " << fmt2str.at(style) << "\n\n";
    makeWriter(tokens, molID, atomtypemap)(file, step);

    if(bonds && !bondlist.empty()){
        file << "\nBonds\n\n";
        for(size_t i=0; i!=bondlist.size(); ++i){
            file << (i+1) << ' '
                 << std::get<2>(bondlist[i]) << ' '
                 << (std::get<0>(bondlist[i])+1) << ' '
                 << (std::get<1>(bondlist[i])+1) << '\n';
        }
    }

    if(angles && !anglelist.empty()){
        file << "\nAngles\n\n";
        for(size_t i=0; i!=anglelist.size(); ++i){
            file << (i+1) << ' '
                 << std::get<3>(anglelist[i]) << ' '
                 << (std::get<0>(anglelist[i])+1) << ' '
                 << (std::get<1>(anglelist[i])+1) << ' '
                 << (std::get<2>(anglelist[i])+1) << '\n';
        }
    }

    if(dihedrals && !dihedlist.empty()){
        file << "\nDihedrals\n\n";
        for(size_t i=0; i!=dihedlist.size(); ++i){
            file << (i+1) << ' '
                 << dihedtypemap.at(dihedtypelist[i]) << ' '
                 << (std::get<0>(dihedlist[i])+1) << ' '
                 << (std::get<1>(dihedlist[i])+1) << ' '
                 << (std::get<2>(dihedlist[i])+1) << ' '
                 << (std::get<3>(dihedlist[i])+1) << '\n';
        }
    }

    if(impropers && !improplist.empty()){
        file << "\nImpropers\n\n";
        for(size_t i=0; i!=improplist.size(); ++i){
            file << (i+1) << ' '
                 << improptypemap.at(improptypelist[i]) << ' '
                 << (std::get<0>(improplist[i])+1) << ' '
                 << (std::get<1>(improplist[i])+1) << ' '
                 << (std::get<2>(improplist[i])+1) << ' '
                 << (std::get<3>(improplist[i])+1) << '\n';
        }
    }
    file << '\n';
    return true;
}

static IO::BasePreset makePreset()
{
    return {&IO::LmpInput,
        {{"style", static_cast<uint>(LMPAtomStyle::Atomic)},
         {"bonds", false},
         {"angles", false},
         {"dihedrals", false},
         {"impropers", false},
    }};
}

const IO::Plugin IO::LmpInput =
{
    "Lammps Data File",
    "lmp",
    "lmp",
    &LmpInpParser,
    &LmpInpWriter,
    nullptr,
    &makePreset
};
