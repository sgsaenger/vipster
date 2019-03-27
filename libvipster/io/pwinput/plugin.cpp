#include "plugin.h"
#include "crystal.h"
#include "io/util.h"

#include "tinyexpr.h"

#include <sstream>
#include <iomanip>
#include <cctype>

using namespace Vipster;

struct CellInp{
    enum CellFmt{None, Alat, Bohr, Angstrom};
    CellFmt fmt{CellFmt::None};
    Mat cell;
};

void parseNamelist(std::string name, std::ifstream& file, IO::PWParam& p)
{
    const std::map<std::string, IO::PWParam::Namelist IO::PWParam::*> nlmap = {
        {"&CONTROL", &IO::PWParam::control},
        {"&SYSTEM", &IO::PWParam::system},
        {"&ELECTRONS", &IO::PWParam::electrons},
        {"&IONS", &IO::PWParam::ions},
        {"&CELL", &IO::PWParam::cell},
    };

    auto nlp = nlmap.find(name);
    if (nlp == nlmap.end()) throw IO::Error("Unknown namelist");
    IO::PWParam::Namelist &nl = p.*(nlp->second);

    std::string line, key;
    size_t beg, end, quote_end;
    const std::string keysep{"= \t\r"};
    const std::string valsep{"=, \t\r"};
    while (std::getline(file, line)) {
        if (line[0] == '/') return;
        if (line[0] == '!') continue;
        end = 0;
        while((beg = line.find_first_not_of(valsep, end)) != line.npos) {
            end = line.find_first_of(keysep, beg);
            key = line.substr(beg, end-beg);
            beg = line.find_first_not_of(valsep, end);
            end = line.find_first_of(valsep, beg);
            quote_end = (end == line.npos) ? line.length()-1 : end-1;
            while ((line[beg] == '"' && line[quote_end] != '"') ||
                   (line[beg] == '\'' && line[quote_end] != '\'')) {
                end = line.find_first_of(valsep, end);
                quote_end = (end == line.npos) ? line.length()-1 : end-1;
            }
            nl[key] = std::string{line, beg, end-beg};
        }
    }
    throw IO::Error("Error in Namelist-parsing");
}

void parseSpecies(std::ifstream& file, Molecule& m, IO::PWParam& p)
{
    auto dataentry = p.system.find("ntyp");
    if (dataentry == p.system.end()) throw IO::Error("ntyp not specified");
    int ntyp = std::stoi(dataentry->second);
    p.system.erase(dataentry);

    std::string line;
    for(int i=0; i<ntyp; ++i) {
        std::getline(file, line);
        while (line[0] == '!' || line[0] == '#') std::getline(file, line);
        std::string name, mass, pwpp;
        std::stringstream linestream{line};
        linestream >> name >> mass >> pwpp;
        if(linestream.fail()) throw IO::Error("Failed to parse species");
        PseEntry &type = (*m.pse)[name];
        type.m = std::stof(mass);
        type.PWPP = pwpp;
    }
}

void parseCoordinates(std::string name, std::ifstream& file,
                      Molecule& m, IO::PWParam& p)
{
    auto dataentry = p.system.find("nat");
    if (dataentry == p.system.end()) throw IO::Error("nat not specified");
    size_t nat = std::stoul(dataentry->second);
    p.system.erase(dataentry);
    Step &s = m.getStep(0);

    size_t pos = name.find_first_of(' ');
    pos = name.find_first_not_of(' ', pos);
    size_t pos2 = name.find_last_not_of(' ');
    if (pos2 != (pos-1)) {
        auto fmt = std::string{name, pos, pos2};
        if(fmt.find("ALAT") != fmt.npos) {
            s.setFmt(AtomFmt::Alat);
        }else if(fmt.find("BOHR") != fmt.npos) {
            s.setFmt(AtomFmt::Bohr);
        }else if(fmt.find("ANGSTROM") != fmt.npos) {
            s.setFmt(AtomFmt::Angstrom);
        }else if(fmt.find("CRYSTAL") != fmt.npos) {
            s.setFmt(AtomFmt::Crystal);
        }else if(fmt.find("CRYSTAL_SG") != fmt.npos) {
            //TODO
            throw IO::Error("CRYSTAL_SG format not yet implemented");
//            s.setFmt(AtomFmt::Crystal);
        }else{
            throw IO::Error("Unknown atom format: "+fmt);
        }
    } else {
        s.setFmt(AtomFmt::Alat);
    }

    s.newAtoms(nat);
    std::string line, coord_expr;
    int err_pos;
    for (auto& at: s) {
        std::getline(file, line);
        while(line[0]=='!' || line[0]=='#'){
             std::getline(file, line);
        }
        std::stringstream linestream{line};
        linestream >> at.name;
        for(size_t i=0; i<3; ++i){
            linestream >> coord_expr;
            at.coord[i] = static_cast<float>(te_interp(coord_expr.c_str(), &err_pos));
            if(err_pos){
                throw IO::Error("Error parsing atom: "+coord_expr);
            }
        }
        if (linestream.fail()) {
            throw IO::Error{"Failed to parse atom"};
        }
        bool x{1},y{1},z{1};
        linestream >> x >> y >> z;
        at.properties->flags[AtomFlag::FixX] = !x;
        at.properties->flags[AtomFlag::FixY] = !y;
        at.properties->flags[AtomFlag::FixZ] = !z;
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
        m.getKPoints().active = KPoints::Fmt::MPG;
        KPoints::MPG &mpg = m.getKPoints().mpg;
        std::stringstream linestream{line};
        linestream >> mpg.x >> mpg.y >> mpg.z >> mpg.sx >> mpg.sy >> mpg.sz;
        if (linestream.fail()) throw IO::Error("Failed to parse automatic K-Points");
    } else {
        m.getKPoints().active = KPoints::Fmt::Discrete;
        if(name.find("CRYSTAL") != name.npos) {
            m.getKPoints().discrete.properties = KPoints::Discrete::crystal;
        }
        if(name.find("_B") != name.npos) {
            m.getKPoints().discrete.properties =
                    static_cast<KPoints::Discrete::Properties>(
                        m.getKPoints().discrete.properties |
                        KPoints::Discrete::band);
        }
        std::string line;
        std::getline(file, line);
        size_t nk;
        std::stringstream{line} >> nk;
        m.getKPoints().discrete.kpoints.resize(nk);
        for(size_t k=0; k!=nk; ++k){
            std::getline(file, line);
            KPoints::Discrete::Point &kp = m.getKPoints().discrete.kpoints[k];
            std::stringstream{line} >> kp.pos[0] >> kp.pos[1] >> kp.pos[2] >> kp.weight;
        }
    }
}

void parseCell(std::string name, std::ifstream& file, CellInp &cell)
{
    std::string line;
    for(size_t i=0; i<3; ++i){
        std::getline(file, line);
        while(line[0]=='!' || line[0]=='#') std::getline(file, line);
        std::stringstream linestream{line};
        linestream >> cell.cell[i][0] >> cell.cell[i][1] >> cell.cell[i][2];
        if (linestream.fail()) throw IO::Error("Failed to parse CELL_PARAMETERS");
    }
    if (name.find("BOHR") != name.npos) cell.fmt = CellInp::Bohr;
    else if (name.find("ANGSTROM") != name.npos) cell.fmt = CellInp::Angstrom;
    else cell.fmt = CellInp::Alat;
}

void createCell(Molecule &m, IO::PWParam &p, CellInp &cell)
{
    Step &s = m.getStep(0);
    IO::PWParam::Namelist& sys = p.system;
    // make sure that relative coordinates are not changed by "scaling" when setting cdm
    auto ibr = sys.find("ibrav");
    if(ibr == sys.end()){
        throw IO::Error("ibrav needs to be specified");
    }
    sys.erase(ibr);
    auto ibrav = std::stoi(ibr->second);
    if(ibrav == 0){
        bool scale = (s.getFmt() >= AtomFmt::Crystal);
        switch (cell.fmt) {
        case CellInp::Bohr:
            s.setCellDim(1, CdmFmt::Bohr, scale);
            break;
        case CellInp::Angstrom:
            s.setCellDim(1, CdmFmt::Angstrom, scale);
            break;
        case CellInp::Alat:
            {
                auto celldm = sys.find("celldm(1)");
                auto cellA = sys.find("A");
                if ((celldm != sys.end()) && (cellA == sys.end())) {
                    s.setCellDim(std::stof(celldm->second), CdmFmt::Bohr, scale);
                    sys.erase(celldm);
                } else if ((celldm == sys.end()) && (cellA != sys.end())) {
                    s.setCellDim(std::stof(cellA->second), CdmFmt::Angstrom, scale);
                    sys.erase(cellA);
                } else if ((celldm == sys.end()) && (cellA == sys.end())) {
                    s.setCellDim(1, CdmFmt::Bohr, scale);
                    break;
                } else {
                    throw IO::Error("Specify either celldm or A,B,C, but not both!");
                }
            }
            break;
        case CellInp::None:
            throw IO::Error("ibrav=0, but no CELL_PARAMETERS were given");
        }
        s.setCellVec(cell.cell, s.getFmt() == AtomFmt::Crystal);
    }else{
        CdmFmt cdmFmt;
        auto celldm = sys.find("celldm(1)");
        auto cellA = sys.find("A");
        if ((celldm != sys.end()) && (cellA == sys.end())) {
            cdmFmt = CdmFmt::Bohr;
            sys.erase(celldm);
        } else if ((celldm == sys.end()) && (cellA != sys.end())) {
            cdmFmt = CdmFmt::Angstrom;
            sys.erase(cellA);
        } else {
            throw IO::Error("Specify either celldm or A,B,C, but not both!");
        }
        Vec axes{}, angles{};
        if(cdmFmt == CdmFmt::Bohr){
            axes[0] = stof(celldm->second);
            auto celldm2 = sys.find("celldm(2)");
            if(celldm2 != sys.end()) axes[1] = stof(celldm2->second);
            auto celldm3 = sys.find("celldm(3)");
            if(celldm3 != sys.end()) axes[2] = stof(celldm3->second);
            auto celldm4 = sys.find("celldm(4)");
            if(celldm4 != sys.end()) angles[0] = stof(celldm4->second);
            auto celldm5 = sys.find("celldm(5)");
            if(celldm5 != sys.end()) angles[1] = stof(celldm5->second);
            auto celldm6 = sys.find("celldm(6)");
            if(celldm6 != sys.end()) angles[2] = stof(celldm6->second);
        }else{
            axes[0] = stof(cellA->second);
            auto cellB = sys.find("B");
            if(cellB != sys.end()) axes[1] = stof(cellB->second)/axes[0];
            auto cellC = sys.find("C");
            if(cellC != sys.end()) axes[2] = stof(cellC->second)/axes[0];
            auto cosAB = sys.find("cosAB");
            if(cosAB != sys.end()) angles[0] = stof(cosAB->second);
            auto cosAC = sys.find("cosAC");
            if(cosAC != sys.end()) angles[1] = stof(cosAC->second);
            auto cosBC = sys.find("cosBC");
            if(cosBC != sys.end()) angles[2] = stof(cosBC->second);
        }
        s.setCellDim(axes[0], cdmFmt, s.getFmt() >= AtomFmt::Crystal);
        auto cell = makeBravais(ibrav, axes, angles);
        s.setCellVec(cell, s.getFmt() == AtomFmt::Crystal);
    }
}

void parseCard(std::string name, std::ifstream& file,
               Molecule& m, IO::PWParam& p,
               CellInp &cell)
{
    if (name.find("ATOMIC_SPECIES") != name.npos) parseSpecies(file, m, p);
    else if (name.find("ATOMIC_POSITIONS") != name.npos) parseCoordinates(name, file, m, p);
    else if (name.find("K_POINTS") != name.npos) parseKPoints(name, file, m);
    else if (name.find("CELL_PARAMETERS") != name.npos) parseCell(name, file, cell);
    else if (name.find("OCCUPATIONS") != name.npos) throw IO::Error("OCCUPATIONS not implemented");
    else if (name.find("CONSTRAINTS") != name.npos) throw IO::Error("CONSTRAINTS not implemented");
    else if (name.find("ATOMIC_FORCES") != name.npos) throw IO::Error("ATOMIC_FORCES not implemented");
}

IO::Data PWInpParser(const std::string& name, std::ifstream &file)
{
    IO::Data d{};
    Molecule &m = d.mol;
    m.setName(name);
    m.newStep();
    d.param = std::make_unique<IO::PWParam>(name);
    IO::PWParam &p = *static_cast<IO::PWParam*>(d.param.get());
    CellInp cell{};

    std::string buf, line;
    while (std::getline(file, buf)) {
        line = IO::trim(buf);
        if (line.empty() || line[0] == '!' || line[0] == '#') continue;
        for (auto &c: line) c = static_cast<char>(std::toupper(c));
        if (line[0] == '&') parseNamelist(line, file, p);
        else parseCard(line, file, m, p, cell);
    }

    createCell(m, p, cell);

    return d;
}

bool PWInpWriter(const Molecule& m, std::ofstream &file,
                 const IO::BaseParam *const p,
                 const IO::BaseConfig *const c,
                 IO::State state)
{
    const auto& s = m.getStep(state.index);
    const auto *pp = dynamic_cast<const IO::PWParam*>(p);
    if(!pp) throw IO::Error("PWI-Writer needs PWScf parameter set");
    const auto *cc = dynamic_cast<const IO::PWConfig*>(c);
    if(!cc) throw IO::Error("PWI-Writer needs PWScf configuration preset");
    std::vector<std::pair<std::string, const IO::PWParam::Namelist*>>
            outNL = {{"control", &pp->control},
                     {"system", &pp->system},
                     {"electrons", &pp->electrons}};
    auto calc = pp->control.find("calculation");
    if(calc != pp->control.end()){
        if(calc->second == "'vc-relax'"){
            outNL.push_back({"ions", &pp->ions});
            outNL.push_back({"cell", &pp->cell});
        }else if(calc->second == "'relax'"){
            outNL.push_back({"ions", &pp->ions});
        }
    }
    for(auto &nl: outNL){
        file << '&' << nl.first << '\n';
        if(nl.second == &pp->system){
            file << " ibrav = 0\n";
            file << " nat = " << s.getNat() << '\n';
            file << " ntyp = " << s.getNtyp() << '\n';
            auto cell_fmt = (cc->cell == IO::PWConfig::CellFmt::Current) ?
                        state.cell_fmt : // use from GUI/CLI
                        static_cast<CdmFmt>(cc->cell); // use from Step
            if(cell_fmt == CdmFmt::Bohr){
                file << " celldm(1) = " << s.getCellDim(cell_fmt) << '\n';
            }else{
                file << " A = " << s.getCellDim(cell_fmt) << '\n';
            }
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
             << std::right << std::setw(9) << e.m << ' ';
        if(e.PWPP.empty()){
            file << t << settings.PWPP.val << '\n';
        }else{
            file << e.PWPP << '\n';
        }
    }
    const std::array<std::string, 4> atfmt = {{"bohr", "angstrom", "crystal", "alat"}};
    auto atom_fmt = (cc->atoms == IO::PWConfig::AtomFmt::Current) ?
                state.atom_fmt : // use from GUI/CLI
                static_cast<AtomFmt>(cc->atoms); // use from Step
    file << "\nATOMIC_POSITIONS " << atfmt[static_cast<size_t>(atom_fmt)] << '\n'
         << std::fixed << std::setprecision(5);
    AtomFlags fixComp{};
    fixComp[AtomFlag::FixX] = true;
    fixComp[AtomFlag::FixY] = true;
    fixComp[AtomFlag::FixZ] = true;
    for (const auto& at: s.asFmt(atom_fmt)) {
        file << std::left << std::setw(3) << at.name
             << ' ' << std::right << std::setw(10) << at.coord[0]
             << ' ' << std::right << std::setw(10) << at.coord[1]
             << ' ' << std::right << std::setw(10) << at.coord[2];
        if((at.properties->flags & fixComp).any()){
            file << ' ' << !at.properties->flags[AtomFlag::FixX]
                 << ' ' << !at.properties->flags[AtomFlag::FixY]
                 << ' ' << !at.properties->flags[AtomFlag::FixZ];
        }
        file << '\n';
    }
    file << "\nK_POINTS " << std::defaultfloat;
    const KPoints& k = m.getKPoints();
    const std::array<std::string, 6> kdprop =
        {{"tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c"}};
    switch(k.active){
    case KPoints::Fmt::Gamma:
        file << "gamma\n";
        break;
    case KPoints::Fmt::MPG:
        file << "automatic\n"
             << k.mpg.x << ' ' << k.mpg.y << ' ' << k.mpg.z << ' '
             << k.mpg.sx << ' ' << k.mpg.sy << ' ' << k.mpg.sz << '\n';
        break;
    case KPoints::Fmt::Discrete:
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

static std::unique_ptr<IO::BaseParam> makeParam(const std::string& name)
{
    return std::make_unique<IO::PWParam>(name);
}

static std::unique_ptr<IO::BaseConfig> makeConfig(const std::string& name)
{
    return std::make_unique<IO::PWConfig>(name);
}

const IO::Plugin IO::PWInput =
{
    "PWScf Input File",
    "pwi",
    "pwi",
    IO::Plugin::Param|IO::Plugin::Config,
    &PWInpParser,
    &PWInpWriter,
    &makeParam,
    &makeConfig
};
