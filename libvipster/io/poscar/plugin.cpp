#include "plugin.h"

#include <sstream>

using namespace Vipster;

IO::Data PoscarParser(const std::string& name, std::ifstream &file){
    IO::Data d{};
    Molecule &m = d.mol;
    m.setName(name);
    Step &s = m.newStep();
    s.setFmt(AtomFmt::Crystal);

    std::string line, tmp;
    size_t n;
    std::stringstream ss;
    std::vector<std::string> types;
    std::vector<size_t> n_per_type;
    bool selective{false};
    float alat;
    auto readline = [&](){
        if(!std::getline(file, line)){
            throw IO::Error{"POSCAR file is missing necessary lines."};
        }
        ss = std::stringstream{line};
    };
    // line 1: comment
    readline();
    s.setComment(line);
    /* line 2: scaling
     * positive: direct factor
     * negative: target volume
     */
    readline();
    ss >> tmp;
    alat = std::stof(tmp);
    // lines 3-5: cell vectors
    Mat cell{};
    for(size_t i=0; i<3; ++i){
        readline();
        ss >> cell[i][0] >> cell[i][1] >> cell[i][2];
        if(ss.fail()){
            throw IO::Error{"Invalid cell vectors in POSCAR file"};
        }
    }
    s.setCellDim(1, CdmFmt::Angstrom);
    s.setCellVec(cell);
    // atom types (optional), then number of atom per type
    readline();
    ss >> tmp;
    if(ss.fail()){
        throw IO::Error{"POSCAR file is missing type information"};
    }
    if(!std::isdigit(tmp[0])){
        // atom types
        types.push_back(tmp);
        while(!ss.eof()){
            ss >> tmp;
            types.push_back(tmp);
        }
        readline();
        ss >> tmp;
    }
    // register first number of atoms that has been prefetched
    n_per_type.push_back(std::stoul(tmp));
    while(!ss.eof()){
        // register rest of cumulative numbers of atoms
        ss >> n;
        n_per_type.push_back(n+n_per_type.back());
    }
    // make sure number of types add up
    if(!types.empty()){
        if(types.size() != n_per_type.size()){
            throw IO::Error{"Mismatching number of atom types in POSCAR file"};
        }
    }else{
        types.resize(n_per_type.size());
        std::iota(types.begin(), types.end(), 1);
    }
    // optional: selective-keyword
    // mandatory: selecting between cartesian or 'direct' aka crystal coordinates
    readline();
    if(std::tolower(line[0]) == 's'){
        selective = true;
        readline();
    }
    if(std::tolower(line[0]) == 'c' || std::tolower(line[0]) == 'k'){
        s.setFmt(AtomFmt::Angstrom);
    }
    // atom coordinates: 3 float columns, if selective 3 columns with T or F
    n = 0;
    std::array<char, 3> sel;
    for(size_t i=0; i<n_per_type.back(); ++i){
        s.newAtom(types[n]);
        if(s.getNat() == n_per_type[n]) n++;
        readline();
        auto at = s[i];
        ss >> at.coord[0] >> at.coord[1] >> at.coord[2];
        if(selective){
            ss >> sel[0] >> sel[1] >> sel[2];
            if(sel[0] == 'F') at.properties->flags[AtomFlag::FixX] = true;
            if(sel[1] == 'F') at.properties->flags[AtomFlag::FixY] = true;
            if(sel[2] == 'F') at.properties->flags[AtomFlag::FixZ] = true;
        }
    }

    // finally, scale coordinates with scaling factor
    if(alat < 0){
        // target volume given, calculate real alat here
        float vol = Mat_det(cell);
        alat = -alat/vol;
    }
    s.setCellDim(alat, CdmFmt::Angstrom, true);

    return d;
}

bool PoscarWriter(const Molecule& m, std::ofstream &file,
                  const IO::BaseParam* const, const IO::BaseConfig* const c,
                  IO::State state)
{
    const auto * const cc = dynamic_cast<const IO::PoscarConfig*const>(c);
    const Step& s = m.getStep(state.index).asFmt(cc->cartesian ?
                                                 AtomFmt::Angstrom : AtomFmt::Crystal);
    file << s.getComment() << '\n';
    file << s.getCellDim(CdmFmt::Angstrom) << '\n';
    file << std::fixed << std::setprecision(10);
    for(size_t i=0; i<3; ++i){
        const auto& v = s.getCellVec()[i];
        file << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
    }
    // collect atoms sorted by their type
    std::map<std::string, std::vector<size_t>> types;
    for(auto it=s.begin(); it!=s.end(); ++it){
        types[it->name].push_back(it.getIdx());
    }
    // print out names of types
    for(const auto& type: types){
        file << ' ' << type.first;
    }
    file << '\n';
    // print out number of atoms per type
    for(const auto& type: types){
        file << ' ' << type.second.size();
    }
    file << '\n';
    // print config options
    if(cc->selective){
        file << "Selective\n";
    }
    if(cc->cartesian){
        file << "Cartesian\n";
    }else{
        file << "Direct\n";
    }
    // print atoms
    auto it = s.begin();
    if(cc->selective){
        for(const auto& type: types){
            for(const auto& idx: type.second){
                it += idx-it.getIdx();
                file << ' ' << it->coord[0]
                     << ' ' << it->coord[1]
                     << ' ' << it->coord[2]
                     << ' ' << (it->properties->flags[AtomFlag::FixX] ? 'F' : 'T')
                     << ' ' << (it->properties->flags[AtomFlag::FixY] ? 'F' : 'T')
                     << ' ' << (it->properties->flags[AtomFlag::FixZ] ? 'F' : 'T')
                     << '\n';
            }
        }
    }else{
        for(const auto& type: types){
            for(const auto& idx: type.second){
                it += idx-it.getIdx();
                file << ' ' << it->coord[0]
                     << ' ' << it->coord[1]
                     << ' ' << it->coord[2] << '\n';
            }
        }
    }
    return true;
}

static std::unique_ptr<IO::BaseConfig> makeConfig(const std::string& name)
{
    return std::make_unique<IO::PoscarConfig>(name);
}

const IO::Plugin IO::Poscar =
{
    "VASP POSCAR",
    "POSCAR",
    "pos",
    IO::Plugin::Config,
    &PoscarParser,
    &PoscarWriter,
    nullptr,
    &makeConfig,
};
