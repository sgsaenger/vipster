#include "ioplugins/pwinput.h"

#include <sstream>
#include <iomanip>

using namespace Vipster;

enum class CellFmt{None, Alat, Bohr, Angstrom};

void parseNamelist(std::string name, std::ifstream& file, IO::PWParam& p)
{
    const std::map<std::string, IO::PWNamelist IO::PWParam::*> nlmap = {
        {"&CONTROL", &IO::PWParam::control},
        {"&SYSTEM", &IO::PWParam::system},
        {"&ELECTRONS", &IO::PWParam::electrons},
        {"&IONS", &IO::PWParam::ions},
        {"&CELL", &IO::PWParam::cell},
    };
    auto nlp = nlmap.find(name);
    if (nlp == nlmap.end()) throw IOError("Unknown namelist");
    IO::PWNamelist &nl = p.*(nlp->second);

    std::string line, key;
    size_t beg, end, quote_end;
    const std::string keysep{"= "};
    const std::string valsep{"=, "};
    while (std::getline(file, line)) {
        if (line[0] == '/') return;
        if (line[0] == '!') continue;
        end = 0;
        while((beg = line.find_first_not_of(keysep, end)) != line.npos) {
            end = line.find_first_of(keysep, beg);
            key = std::string{line, beg, end-1};
            beg = line.find_first_not_of(valsep, end);
            end = line.find_first_of(valsep, beg);
            quote_end = (end == line.npos) ? line.length()-1 : end;
            while ((line[beg] == '"' && line[quote_end] != '"') ||
                   (line[beg] == '\'' && line[quote_end] != '\'')) {
                line = std::string{line, end};
                end = line.find_first_of(valsep, end);
            }
            nl[key] = std::string{line, beg, end};
        }
    }
    throw IOError("Error in Namelist-parsing");
}

void parseSpecies(std::ifstream& file, IO::PWData& d)
{
    auto dataentry = d.data.system.find("ntyp");
    if (dataentry == d.data.system.end()) throw IOError("ntyp not specified");
    int ntyp = std::stoi(dataentry->second);
    d.data.system.erase(dataentry);

    std::string line;
    for(int i=0; i<ntyp; ++i) {
        std::getline(file, line);
        while (line[0] == '!' || line[0] == '#') std::getline(file, line);
        std::string name, mass, pwpp;
        std::stringstream linestream{line};
        linestream >> name >> mass >> pwpp;
        if(linestream.fail()) throw IOError("Failed to parse species");
        PseEntry &type = (*d.mol.pse)[name];
        type.m = std::stof(mass);
        type.PWPP = pwpp;
    }
}

void parseCoordinates(std::string name, std::ifstream& file, IO::PWData& d)
{
    auto dataentry = d.data.system.find("nat");
    if (dataentry == d.data.system.end()) throw IOError("nat not specified");
    size_t nat = std::stoul(dataentry->second);
    d.data.system.erase(dataentry);
    StepProper &s = d.mol.getStep(0);

    const std::map<std::string, AtomFmt> fmtmap = {
        {"ALAT", AtomFmt::Alat},
        {"BOHR", AtomFmt::Bohr},
        {"CRYSTAL", AtomFmt::Crystal},
        {"ANGSTROM", AtomFmt::Angstrom},
        {"CRYSTAL_SG", AtomFmt::Crystal}
    };
    size_t pos = name.find_first_of(' ');
    pos = name.find_first_not_of(' ', pos);
    size_t pos2 = name.find_last_not_of(' ');
    if (pos2 != (pos-1)) {
        auto fmt = std::string{name, pos, pos2};
        if (fmt == "CRYSTAL_SG") throw IOError("CRYSTAL_SG format not implemented");
        auto atfmt = fmtmap.find(fmt);
        if (atfmt == fmtmap.end()) throw IOError("Unknown atom format");
        s.setFmt(atfmt->second);
    } else {
        s.setFmt(AtomFmt::Alat);
    }

    s.newAtoms(nat);
    std::string line;
    for (size_t i=0; i<nat; ++i) {
        std::getline(file, line);
        while(line[0]=='!' || line[0]=='#') std::getline(file, line);
        auto at = s[i];
        std::stringstream linestream{line};
        linestream >> at.name >> at.coord[0] >> at.coord[1] >> at.coord[2];
        if (linestream.fail()) throw IOError{"Failed to parse atom"};
        linestream >> at.fix[0] >> at.fix[1] >> at.fix[2];
    }
}

void parseKPoints(std::string name, std::ifstream& file, Molecule& m)
{
    if (name.find("GAMMA") != name.npos) {
        return;
    } else if (name.find("AUTOMATIC") != name.npos) {
        std::string line;
        std::getline(file, line);
        while(line[0]=='!' || line[0]=='#') std::getline(file, line);
        KPoints::MPG &mpg = m.getKPoints().mpg;
        std::stringstream linestream{line};
        linestream >> mpg.x >> mpg.y >> mpg.z >> mpg.sx >> mpg.sy >> mpg.sz;
        if (linestream.fail()) throw IOError("Failed to parse automatic K-Points");
    } else {
        throw IOError("Discrete K-Points not implemented");
    }
}

void parseCell(std::string name, std::ifstream& file, IO::PWData& d, CellFmt &cellFmt)
{
    auto ibrav = d.data.system.find("ibrav");
    if (ibrav == d.data.system.end()) throw IOError{"ibrav not specified"};
    if(!std::stoi(ibrav->second)) {
        std::string line;
        Mat cell;
        for(int i=0; i<3; ++i){
            std::getline(file, line);
            while(line[0]=='!' || line[0]=='#') std::getline(file, line);
            std::stringstream linestream{line};
            linestream >> cell[i][0] >> cell[i][1] >> cell[i][2];
            if (linestream.fail()) throw IOError("Failed to parse CELL_PARAMETERS");
        }
        Step& step = d.mol.getStep(0);
        step.setCellVec(cell, (step.getFmt()==AtomFmt::Crystal));
        if (name.find("BOHR") != name.npos) cellFmt = CellFmt::Bohr;
        else if (name.find("ANGSTROM") != name.npos) cellFmt = CellFmt::Angstrom;
        else cellFmt = CellFmt::Alat;
    }
}

void createCell(IO::PWData &d, CellFmt &cellFmt)
{
    StepProper &s = d.mol.getStep(0);
    CdmFmt cdmFmt;
    const IO::PWNamelist& sys = d.data.system;
    auto celldm = sys.find("celldm(1)");
    auto cellA = sys.find("A");
    if ((celldm != sys.end()) && (cellA == sys.end())) {
        cdmFmt = CdmFmt::Bohr;
    } else if ((celldm == sys.end()) && (cellA != sys.end())) {
        cdmFmt = CdmFmt::Angstrom;
    } else {
        throw IOError("Specify either celldm or A,B,C, but not both!");
    }
    bool scale = (s.getFmt() >= AtomFmt::Crystal);
    switch (cellFmt) {
    case CellFmt::Bohr:
        s.setCellDim(1, CdmFmt::Bohr, scale);
        break;
    case CellFmt::Angstrom:
        s.setCellDim(1, CdmFmt::Angstrom, scale);
        break;
    case CellFmt::Alat:
        switch (cdmFmt) {
        case CdmFmt::Angstrom:
            s.setCellDim(std::stof(cellA->second), CdmFmt::Angstrom, scale);
            break;
        case CdmFmt::Bohr:
            s.setCellDim(std::stof(celldm->second), CdmFmt::Bohr, scale);
            break;
        }
        break;
    case CellFmt::None:
        auto ibrav = d.data.system.find("ibrav");
        if (ibrav == d.data.system.end()) throw IOError{"ibrav not specified"};
        if(!std::stoi(ibrav->second)) throw IOError("ibrav=0, but no CELL_PARAMETERS were given");
        //TODO
        throw IOError("Creating Cells based on ibrav not supported yet");
        break;
    }
}

void parseCard(std::string name, std::ifstream& file, IO::PWData& d, CellFmt &cellFmt)
{
    if (name.find("ATOMIC_SPECIES") != name.npos) parseSpecies(file, d);
    else if (name.find("ATOMIC_POSITIONS") != name.npos) parseCoordinates(name, file, d);
    else if (name.find("K_POINTS") != name.npos) parseKPoints(name, file, d.mol);
    else if (name.find("CELL_PARAMETERS") != name.npos) parseCell(name, file, d, cellFmt);
    else if (name.find("OCCUPATIONS") != name.npos) throw IOError("OCCUPATIONS not implemented");
    else if (name.find("CONSTRAINTS") != name.npos) throw IOError("CONSTRAINTS not implemented");
    else if (name.find("ATOMIC_FORCES") != name.npos) throw IOError("ATOMIC_FORCES not implemented");
}

std::shared_ptr<IO::BaseData> PWInpParser(std::string name, std::ifstream &file)
{
    auto d = std::make_shared<IO::PWData>();
    Molecule &m = d->mol;
    m.setName(name);
    m.newStep();
    IO::PWParam &p = d->data;
    CellFmt cellFmt = CellFmt::None;

    std::string line;
    while (std::getline(file, line)) {
        if (!line[0] || line[0] == ' ' || line[0] == '!' || line[0] == '#') continue;
        for (auto &c: line) c = std::toupper(c);
        if (line[0] == '&') parseNamelist(line, file, p);
        else parseCard(line, file, *d, cellFmt);
    }

    createCell(*d, cellFmt);

    return d;
}

bool PWInpWriter(const Molecule& m, std::ofstream &file, const IO::BaseParam* p)
{
    auto& s = m.getStep(0);
    auto *pp = dynamic_cast<const IO::PWParam*>(p);
    if(!pp) throw IOError("PWI-Writer needs PWScf parameter set");
    std::vector<std::pair<std::string, const std::map<std::string,std::string>*>>
            outNL = {{"control", &pp->control},
                     {"system", &pp->system},
                     {"electrons", &pp->electrons}};
    auto c = pp->control.find("calculation");
    if(c != pp->control.end()){
        if(c->second == "'vc-relax'"){
            outNL.push_back({"ions", &pp->ions});
            outNL.push_back({"cell", &pp->cell});
        }else if(c->second == "'relax'"){
            outNL.push_back({"ions", &pp->ions});
        }
    }
    for(auto &nl: outNL){
        file << '&' << nl.first << '\n';
        if(nl.second == &pp->system){
            file << " nat = " << s.getNat() << '\n';
            file << " ntyp = " << s.getNtyp() << '\n';
            //TODO: decide for A when ...?
            file << " celldm(1) = " << s.getCellDim(CdmFmt::Bohr) << '\n';
        }
        for(auto& e: *nl.second){
            file << ' ' << e.first << " = " << e.second << '\n';
        }
        file << "/\n\n";
    }
    file << "ATOMIC_SPECIES\n"
         << std::fixed << std::setprecision(5);
    for(auto &t: s.getTypes()){
        auto e = (*s.pse)[t];
        file << std::left << std::setw(3) << t << ' '
             << std::right << std::setw(9) << e.m << ' '
             << e.PWPP << '\n';
    }
    const std::array<std::string, 4> atfmt = {{"bohr", "angstrom", "crystal", "alat"}};
    file << "\nATOMIC_POSITION " << atfmt[(int)s.getFmt()] << '\n'
         << std::fixed << std::setprecision(5);
    for (const Atom& at: s) {
        file << std::left << std::setw(3) << at.name << ' '
             << std::right << std::setw(10) << at.coord[0] << ' '
             << std::right << std::setw(10) << at.coord[1] << ' '
             << std::right << std::setw(10) << at.coord[2] << '\n';
    }
    file << "\nK_POINTS " << std::defaultfloat;
    const KPoints& k = m.getKPoints();
    const std::array<std::string, 6> kdprop =
        {{"tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c"}};
    switch(k.active){
    case KPointFmt::Gamma:
        file << "gamma\n";
        break;
    case KPointFmt::MPG:
        file << "automatic\n"
             << k.mpg.x << ' ' << k.mpg.y << ' ' << k.mpg.z << ' '
             << k.mpg.sx << ' ' << k.mpg.sy << ' ' << k.mpg.sz << '\n';
        break;
    case KPointFmt::Discrete:
        file << kdprop[k.discrete.properties] << '\n'
             << k.discrete.kpoints.size() << '\n';
        for(auto &kd: k.discrete.kpoints){
            file << kd.pos[0] << ' ' << kd.pos[1] << ' '
                 << kd.pos[2] << ' ' << kd.weight << '\n';
        }
    }
    file << "\nCELL_PARAMETERS\n" << std::fixed << std::setprecision(5);
    for(auto &v: s.getCellVec()){
        file << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
    }
    return true;
}

const IOPlugin IO::PWInput =
{
    "PWScf Input File",
    "pwi",
    "pwi",
    &PWInpParser,
    &PWInpWriter
};
