#include "plugin.h"
#include <sstream>

using namespace Vipster;

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
    {"charge", {{lmpTok::charge, lmpTok::type, lmpTok::pos}}},
    {"dipole", {{lmpTok::charge, lmpTok::type, lmpTok::pos}}},
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

std::vector<lmpTok> getFmtGuess(std::ifstream& file, size_t nat){
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
                auto f = stof(at[col]);
                if(f != floorf(f)){
                    return false;
                }
            }
            return true;
        }catch(...){
            return false;
        }
    };
    auto checkDummy = [&atoms](size_t col){
        std::set<float> s, n{0.f};
        for (auto& at: atoms){
            s.insert(stof(at[col]));
        }
        return (s.size()==1)&&(s==n);
    };
    if(narg == 5){
        // only one possible setup (atomic)
        return {lmpTok::type, lmpTok::pos};
    }
    if(narg == 6){
        if(checkInt(2) && checkDummy(2)){
            // angle/molecular have molID, then type
            return {lmpTok::ignore, lmpTok::type, lmpTok::pos};
        }
        // parse as charge, even though this is most likely wrong
        return {lmpTok::type, lmpTok::charge, lmpTok::pos};
    }
    std::vector<lmpTok> parser{};
    size_t col{1}, poscoord{narg-3};
    /* assume:
     * - trailing int-columns are image-flags
     * - three cols before image-flags are position
     * - second or third col are atomtype
     * - first col between type and pos is charge (if present)
     */
    if (checkInt(2) && !checkDummy(2)) {
        // assume col1 is molID, probably fails for ellipsoid
        parser.push_back(lmpTok::ignore);
        col++;
    }
    parser.push_back(lmpTok::type);
    for(size_t img=narg-1; img>=std::min(narg-3,static_cast<size_t>(5)); ++img){
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
    return [fmt](std::ifstream& file, Step& s){
        std::string line{}, dummy{};
        for (auto& at:s) {
            std::getline(file, line);
            std::stringstream ss{line};
            ss >> dummy;
            for (lmpTok tok: fmt) {
                switch(tok){
                case lmpTok::type:
                    ss >> at.name;
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

std::vector<std::set<size_t>> groupSets(const std::vector<std::set<size_t>>& input){
    std::vector<std::set<size_t>> output{input[0]};
    for(auto it=input.begin()+1; it!=input.end(); ++it){
        bool matched{false};
        for(auto& set:output){
            std::set<size_t> test;
            std::set_intersection(set.begin(), set.end(),
                                  it->begin(), it->end(),
                                  std::inserter(test, test.begin()));
            if(!test.empty()){
                matched = true;
                test.clear();
                std::set_union(set.begin(), set.end(),
                               it->begin(), it->end(),
                               std::inserter(test, test.begin()));
                std::swap(set, test);
            }
        }
        if(!matched){
            output.push_back(*it);
        }
    }
    if(input.size() != output.size()){
        return groupSets(output);
    }else{
        return output;
    }
}

auto makeWriter(const std::vector<lmpTok>& fmt,
                const std::vector<size_t>& molID,
                const std::map<std::string, size_t>& atomtypemap){
    return [&](std::ostream& file, const Step& s){
        for(Step::constIterator it=s.begin(); it!=s.end(); ++it){
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
                         << std::setw(10) << (*it).coord[0]
                         << ' ' << std::setw(10) << (*it).coord[1]
                         << ' ' << std::setw(10) << (*it).coord[2];
                    break;
                case lmpTok::ignore:
                    break;
                }
            }
            file << '\n' << std::left;
        }
    };
}

IO::Data LmpInpParser(const std::string& name, std::ifstream &file)
{
    enum class ParseMode{Header,Atoms,Types};

    IO::Data data{};
    data.fmt = IOFmt::LMP;
    Molecule& m = data.mol;
    m.setName(name);
    StepProper& s = m.newStep();
    s.setFmt(AtomFmt::Angstrom);
    s.setCellDim(1, CdmFmt::Angstrom);

    std::string line;
    size_t nat{}, ntype{};
    float t1, t2;
    Mat cell{};
    std::map<std::string, std::string> types{};
    bool hasNames{false};
    while (std::getline(file, line)) {
        if (line.find("atoms") != std::string::npos) {
            std::stringstream ss{line};
            ss >> nat;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to"
                              "parse number of atoms");
            }
            s.newAtoms(nat);
        } else if (line.find("atom types") != std::string::npos) {
            std::stringstream ss{line};
            ss >> ntype;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to"
                              "parse number of types");
            }
        } else if (line.find("xlo xhi") != std::string::npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to"
                              "parse cell X dimension");
            }
            cell[0][0] = t2 - t1;
        } else if (line.find("ylo yhi") != std::string::npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to"
                              "parse cell Y dimension");
            }
            cell[1][1] = t2 - t1;
        } else if (line.find("zlo zhi") != std::string::npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to"
                              "parse cell Z dimension");
            }
            cell[2][2] = t2 - t1;
        } else if (line.find("xy xz yz") != std::string::npos) {
            std::stringstream ss{line};
            ss >> cell[1][0] >> cell[2][0] >> cell[2][1];
            if (ss.fail()) {
                throw IO::Error("Lammps Input: failed to"
                              "parse cell tilt factors");
            }
        } else if (line.find("Masses") != std::string::npos) {
            std::getline(file, line);
            std::string id, name;
            for (size_t i=0; i<ntype; ++i) {
                std::getline(file, line);
                std::stringstream ss{line};
                ss >> id >> t1;
                std::size_t cpos = line.find('#');
                if(cpos != std::string::npos) {
                    // if there's a comment, extract the typename from it
                    hasNames = true;
                    std::stringstream ss{line.substr(cpos+1)};
                    ss >> name;
                    types[id] = name;
                } else {
                    // else just number the types accordingly
                    types[id] = id;
                }
                if (ss.fail()) {
                    throw IO::Error("Lammps Input: failed to parse atom type");
                }
                (*s.pse)[name].m = t1;
            }
        } else if (line.find("Atoms") != std::string::npos) {
            std::vector<lmpTok> fmt{};
            // lookup fixed parser if format is given
            std::size_t cpos = line.find('#');
            if (cpos != std::string::npos) {
                std::string f{};
                std::stringstream{line.substr(cpos+1)} >> f;
                fmt = fmtmap.at(f);
            }
            std::getline(file, line);
            // if no format was given, try to determine a suitable parser
            if (fmt.empty()) {
                fmt = getFmtGuess(file, nat);
            }
            // do the parsing
            makeParser(fmt)(file, s);
        }
    }
    s.setCellVec(cell);
    if (hasNames) {
        for (auto& at: s) {
            at.name = types[at.name];
        }
    }
    return data;
}

bool LmpInpWriter(const Molecule& m, std::ofstream &file,
                  const BaseParam *const,
                  const BaseConfig *const c,
                  IO::State state)
{
    const auto& step = m.getStep(state.index);
    const auto& types = step.getNames();
    const auto *cc = dynamic_cast<const IO::LmpConfig*>(c);
    if(!cc) throw IO::Error("Lammps-Writer needs configuration preset");
    const auto tokens = fmtmap.at(IO::LmpConfig::fmt2str.at(cc->style));
    bool needsMolID = std::find(tokens.begin(), tokens.end(), lmpTok::mol) != tokens.end();

    file << std::fixed << std::setprecision(5);

    // prepare bonds
    std::vector<std::tuple<size_t, size_t>> bondlist;
    std::vector<std::tuple<std::string, std::string>> bondtypelist;
    std::map<std::tuple<std::string, std::string>, size_t> bondtypemap;
    if(cc->bonds || cc->angles || cc->dihedrals || cc->impropers || needsMolID){
        for(const auto& bond: step.getBonds()){
            bondlist.emplace_back(std::min(bond.at1, bond.at2),
                                  std::max(bond.at1, bond.at2));
        }
        std::sort(bondlist.begin(), bondlist.end());
        for(const auto& bond: bondlist){
            bool begin0 = types[std::get<0>(bond)] < types[std::get<1>(bond)];
            if(begin0){
                bondtypelist.emplace_back(types[std::get<0>(bond)], types[std::get<1>(bond)]);
            }else{
                bondtypelist.emplace_back(types[std::get<1>(bond)], types[std::get<0>(bond)]);
            }
            bondtypemap.emplace(bondtypelist.back(), bondtypemap.size()+1);
        }
    }

    // prepare angles
    std::vector<std::tuple<size_t, size_t, size_t>> anglelist;
    std::vector<std::tuple<std::string, std::string, std::string>> angletypelist;
    std::map<std::tuple<std::string, std::string, std::string>, size_t> angletypemap;
    if(cc->angles || cc->dihedrals || cc->impropers){
        for(auto it=bondlist.begin(); it!=bondlist.end(); ++it){
            for(auto it2=it+1; it2!=bondlist.end(); ++it2){
                const auto& a0 = std::get<0>(*it);
                const auto& a1 = std::get<1>(*it);
                const auto& b0 = std::get<0>(*it2);
                const auto& b1 = std::get<1>(*it2);
                if(a0 == b0){
                    anglelist.emplace_back(a1, a0, b1);
                }else if(a1 == b0){
                    anglelist.emplace_back(a0, a1, b1);
                }else if(a1 == b1){
                    anglelist.emplace_back(a0, a1, b0);
                }
                // a0 == b1 impossible because of sorting
            }
        }
        for(const auto& angle: anglelist){
            bool begin0 = types[std::get<0>(angle)] < types[std::get<2>(angle)];
            if(begin0){
                angletypelist.emplace_back(types[std::get<0>(angle)],
                                           types[std::get<1>(angle)],
                                           types[std::get<2>(angle)]);
            }else{
                angletypelist.emplace_back(types[std::get<2>(angle)],
                                           types[std::get<1>(angle)],
                                           types[std::get<0>(angle)]);
            }
            angletypemap.emplace(angletypelist.back(), angletypemap.size()+1);
        }
    }

    // prepare dihedrals
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> dihedlist;
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> dihedtypelist;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, size_t> dihedtypemap;
    if(cc->dihedrals){
        for(auto it=anglelist.begin(); it!=anglelist.end(); ++it){
            for(auto it2=it+1; it2!=anglelist.end(); ++it2){
                const auto& a0 = std::get<0>(*it);
                const auto& a1 = std::get<1>(*it);
                const auto& a2 = std::get<2>(*it);
                const auto& b0 = std::get<0>(*it2);
                const auto& b1 = std::get<1>(*it2);
                const auto& b2 = std::get<2>(*it2);
                if(a1 == b0){
                    if(b1 == a0){
                        dihedlist.emplace_back(a2, a1, a0, b2);
                    }else if(b1 == a2){
                        dihedlist.emplace_back(a0, a1, a2, b2);
                    }
                }else if(a1 == b2){
                    if(b1 == a0){
                        dihedlist.emplace_back(a2, a1, a0, b0);
                    }else if(b1 == a2){
                        dihedlist.emplace_back(a0, a1, a2, b0);
                    }
                }
            }
        }
        for(const auto& dihed: dihedlist){
            const auto& d0 = std::get<0>(dihed);
            const auto& d1 = std::get<1>(dihed);
            const auto& d2 = std::get<2>(dihed);
            const auto& d3 = std::get<3>(dihed);
            bool begin0 = (types[d0] < types[d3]) ||
                          ((types[d0] == types[d3]) && (types[d1] < types[d2]));
            if(begin0){
                dihedtypelist.emplace_back(types[d0], types[d1], types[d2], types[d3]);
            }else{
                dihedtypelist.emplace_back(types[d3], types[d2], types[d1], types[d0]);
            }
            dihedtypemap.emplace(dihedtypelist.back(), dihedtypemap.size()+1);
        }
    }

    // prepare impropers
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> improplist;
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> improptypelist;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, size_t> improptypemap;
    if(cc->impropers){
        for(auto it=anglelist.begin(); it!=anglelist.end(); ++it){
            for(auto it2=it+1; it2!=anglelist.end(); ++it2){
                const auto& a0 = std::get<0>(*it);
                const auto& a1 = std::get<1>(*it);
                const auto& a2 = std::get<2>(*it);
                const auto& b0 = std::get<0>(*it2);
                const auto& b1 = std::get<1>(*it2);
                const auto& b2 = std::get<2>(*it2);
                if(a1 == b1){
                    if(a0 == b0){
                        improplist.emplace_back(a1, a0, a2, b2);
                    }else if(a0 == b2){
                        improplist.emplace_back(a1, a0, a2, b0);
                    }else if(a2 == b0){
                        improplist.emplace_back(a1, a0, a2, b2);
                    }else if(a2 == b2){
                        improplist.emplace_back(a1, a0, a2, b0);
                    }
                }
            }
        }
        for(const auto& improp: improplist){
            const auto& i0 = std::get<0>(improp);
            auto i1 = std::get<1>(improp);
            auto i2 = std::get<2>(improp);
            auto i3 = std::get<3>(improp);
            if(types[i1]>types[i2]){
                std::swap(i1, i2);
            }
            if(types[i1]>types[i3]){
                std::swap(i1,i3);
            }
            if(types[i2]>types[i3]){
                std::swap(i2,i3);
            }
            improptypelist.emplace_back(types[i0], types[i1], types[i2], types[i3]);
            improptypemap.emplace(improptypelist.back(), improptypemap.size()+1);
        }
    }

    // prepare Molecule-IDs
    std::vector<size_t> molID{};
    std::vector<std::set<size_t>> molSets{};
    if(needsMolID){
        molID.resize(step.getNat());
        for(const auto& bond: step.getBonds()){
            molSets.push_back(std::set<size_t>{bond.at1, bond.at2});
        }
        molSets = groupSets(molSets);
        for(size_t i=0; i<molSets.size(); ++i){
            for(const auto& at: molSets[i]){
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
    if(cc->bonds && !bondlist.empty()){
        file << bondlist.size() << " bonds\n"
             << bondtypemap.size() << " bond types\n";
        for(const auto& pair: bondtypemap){
            file << '#' << pair.second << ' '
                 << std::get<0>(pair.first) << ' '
                 << std::get<1>(pair.first) << '\n';
        }
    }
    if(cc->angles && !anglelist.empty()){
        file << anglelist.size() << " angles\n"
             << angletypemap.size() << " angle types\n";
        for(const auto& pair: angletypemap){
            file << '#' << pair.second << ' '
                 << std::get<0>(pair.first) << ' '
                 << std::get<1>(pair.first) << ' '
                 << std::get<2>(pair.first) << '\n';
        }
    }
    if(cc->dihedrals && !dihedlist.empty()){
        file << dihedlist.size() << " dihedrals\n"
             << dihedtypemap.size() << " dihedral types\n";
        for(const auto& pair: dihedtypemap){
            file << '#' << pair.second << ' '
                 << std::get<0>(pair.first) << ' '
                 << std::get<1>(pair.first) << ' '
                 << std::get<2>(pair.first) << ' '
                 << std::get<3>(pair.first) << '\n';
        }
    }
    if(cc->impropers && !improplist.empty()){
        file << improplist.size() << " impropers\n"
             << improptypemap.size() << " improper types\n";
        for(const auto& pair: improptypemap){
            file << '#' << pair.second << ' '
                 << std::get<0>(pair.first) << ' '
                 << std::get<1>(pair.first) << ' '
                 << std::get<2>(pair.first) << ' '
                 << std::get<3>(pair.first) << '\n';
        }
    }

    auto vec = step.getCellVec() * step.getCellDim(CdmFmt::Angstrom);
    if(!float_comp(vec[0][1], 0.f) || !float_comp(vec[0][2], 0.f) || !float_comp(vec[1][2], 0.f)){
        throw IO::Error("Cell vectors must form diagonal or lower triangular matrix for Lammps");
    }
    file << "\n0.0 "
         << vec[0][0] << " xlo xhi\n0.0 "
         << vec[1][1] << " ylo yhi\n0.0 "
         << vec[2][2] << " zlo zhi\n";
    if(!float_comp(vec[1][0], 0.f) || !float_comp(vec[2][0], 0.f) || !float_comp(vec[2][1], 0.f)){
        file << vec[1][0] << ' ' << vec[2][0] << ' ' << vec[2][1] << " xy xz yz\n";
    }

    file << "\nMasses\n\n";
    std::map<std::string, size_t> atomtypemap;
    for(const auto& i: step.getTypes()){
        atomtypemap.emplace(i, atomtypemap.size()+1);
        file << step.pse->at(i).m << " # " << i << '\n';
    }

    file << "\nAtoms # " << IO::LmpConfig::fmt2str.at(cc->style) << "\n\n";
    makeWriter(tokens, molID, atomtypemap)(file, step);

    if(cc->bonds && !bondlist.empty()){
        file << "\nBonds\n\n";
        for(size_t i=0; i!=bondlist.size(); ++i){
            file << (i+1) << ' '
                 << bondtypemap.at(bondtypelist[i]) << ' '
                 << (std::get<0>(bondlist[i])+1) << ' '
                 << (std::get<1>(bondlist[i])+1) << '\n';
        }
    }

    if(cc->angles && !anglelist.empty()){
        file << "\nAngles\n\n";
        for(size_t i=0; i!=anglelist.size(); ++i){
            file << (i+1) << ' '
                 << angletypemap.at(angletypelist[i]) << ' '
                 << (std::get<0>(anglelist[i])+1) << ' '
                 << (std::get<1>(anglelist[i])+1) << ' '
                 << (std::get<2>(anglelist[i])+1) << '\n';
        }
    }

    if(cc->dihedrals && !dihedlist.empty()){
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

    if(cc->impropers && !improplist.empty()){
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

static std::unique_ptr<BaseConfig> makeConfig(const std::string& name)
{
    return std::make_unique<IO::LmpConfig>(name);
}

const IO::Plugin IO::LmpInput =
{
    "Lammps Data File",
    "lmp",
    "lmp",
    IO::Plugin::Config,
    &LmpInpParser,
    &LmpInpWriter,
    nullptr,
    &makeConfig
};
