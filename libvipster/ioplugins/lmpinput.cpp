#include "lmpinput.h"

#include <sstream>

using namespace Vipster;

enum class lmpTok{
    type,
    pos,
    charge,
    ignore
};

const static std::map<std::string, std::vector<lmpTok>> fmtmap{
    {"angle", {{lmpTok::ignore, lmpTok::type, lmpTok::pos}}},
    {"atomic", {{lmpTok::type, lmpTok::pos}}},
    {"body", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"bond", {{lmpTok::ignore, lmpTok::type, lmpTok::pos}}},
    {"charge", {{lmpTok::charge, lmpTok::type, lmpTok::pos}}},
    {"dipole", {{lmpTok::charge, lmpTok::type, lmpTok::pos}}},
    {"dpd", {{lmpTok::type, lmpTok::ignore, lmpTok::pos}}},
    {"edpd", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"mdpd", {{lmpTok::type, lmpTok::pos}}},
    {"tdpd", {{lmpTok::type, lmpTok::pos}}},
    {"electron", {{lmpTok::type, lmpTok::charge, lmpTok::ignore,
                   lmpTok::ignore, lmpTok::pos}}},
    {"ellipsoid", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"full", {{lmpTok::ignore, lmpTok::type, lmpTok::charge, lmpTok::pos}}},
    {"line", {{lmpTok::ignore, lmpTok::type, lmpTok::ignore,
               lmpTok::ignore, lmpTok::pos}}},
    {"meso", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore,
               lmpTok::ignore, lmpTok::pos}}},
    {"molecular", {{lmpTok::ignore, lmpTok::type, lmpTok::pos}}},
    {"peri", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"smd", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::ignore,
              lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"sphere", {{lmpTok::type, lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"template", {{lmpTok::ignore, lmpTok::ignore, lmpTok::ignore,
                   lmpTok::type, lmpTok::pos}}},
    {"tri", {{lmpTok::ignore, lmpTok::type, lmpTok::ignore,
              lmpTok::ignore, lmpTok::pos}}},
    {"wavepacket", {{lmpTok::type, lmpTok::charge, lmpTok::ignore,lmpTok::ignore,
                     lmpTok::ignore, lmpTok::ignore, lmpTok::pos}}},
    {"hybrid", {{lmpTok::type, lmpTok::pos}}}
};

std::vector<lmpTok> getFmtGuess(std::ifstream& file, size_t nat){
    // WILL fail if fmt == tdpd, hybrid, template
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
        std::set<float> s, n{0};
        for (auto& at: atoms){
            s.insert(stof(at[col]));
        }
        return (s.size()==1)&&(s==n);
    };
    if(narg == 5){
        // only one possible setup (atomic, mdpd)
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
                    ss >> at.charge;
                    break;
                case lmpTok::pos:
                    ss >> at.coord[0] >> at.coord[1] >> at.coord[2];
                    break;
                case lmpTok::ignore:
                    ss >> dummy;
                    break;
                }
            }
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

bool LmpInpWriter(const Molecule&, std::ofstream &,
                  const BaseParam *const,
                  const BaseConfig *const c)
{
    auto *lp = dynamic_cast<const IO::LmpConfig*>(c);
    if(lp == nullptr){
        throw IO::Error("Lammps-Writer needs parameter set");
    }
    return false;
}

const IO::Plugin IO::LmpInput =
{
    "Lammps Data File",
    "lmp",
    "lmp",
    IO::Plugin::Config,
    &LmpInpParser,
    nullptr
};

std::unique_ptr<BaseConfig> IO::LmpConfig::copy()
{
    return std::make_unique<IO::LmpConfig>(*this);
}
