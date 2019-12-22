#include "pwinput.h"
#include "crystal.h"
#include "io/util.h"

#include "tinyexpr.h"

#include <sstream>
#include <iomanip>
#include <cctype>

using namespace Vipster;

using NameList = std::map<std::string, std::string>;

static IO::Parameter makeParam()
{
    return {&IO::PWInput, {
            {"&CONTROL", NameList{}},
            {"&SYSTEM", NameList{}},
            {"&ELECTRONS", NameList{}},
            {"&IONS", NameList{}},
            {"&CELL", NameList{}},
            {"PPPrefix", std::string{}},
            {"PPSuffix", std::string{}},
        }};
}

static IO::Preset makePreset()
{
    return {&IO::PWInput,
        {{"atoms", {NamedEnum{4, {"Bohr", "Angstrom", "Crystal", "Alat", "Active"}},
                    "Active: Use the current Step's active atom format\n"
                    "Else: Enforce the selected format"}},
         {"cell", {NamedEnum{2, {"Bohr", "Angstrom", "Active"}},
                    "Active: Match current Step's active atom format (Å if Å, Bohr otherwise)\n"
                    "Else: Enforce the selected format"}}}};
}

struct CellInp{
    enum CellFmt{None, Alat, Bohr, Angstrom};
    CellFmt fmt{CellFmt::None};
    Mat cell;
};

void parseNamelist(std::string name, std::istream& file, IO::Parameter& p)
{
    auto pos = p.find(name);
    if (pos == p.end()) throw IO::Error{"Unknown namelist"};
    auto &nl = std::get<NameList>(pos->second);

    std::string line, key;
    size_t beg, end, quote_end;
    const std::string keysep{"= \t\r"};
    const std::string valsep{"=, \t\r"};
    while (std::getline(file, line)) {
        line = IO::trim(line);
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

void parseSpecies(std::istream& file, Molecule& m, IO::Parameter& p)
{
    auto& sys = std::get<NameList>(p.at("&SYSTEM"));
    auto dataentry = sys.find("ntyp");
    if (dataentry == sys.end()) throw IO::Error("ntyp not specified");
    int ntyp = std::stoi(dataentry->second);
    sys.erase(dataentry);

    std::string line;
    for(int i=0; i<ntyp; ++i) {
        std::getline(file, line);
        line = IO::trim(line);
        while (line[0] == '!' || line[0] == '#'){
            std::getline(file, line);
            line = IO::trim(line);
        }
        std::string name, mass, pwpp;
        std::stringstream linestream{line};
        linestream >> name >> mass >> pwpp;
        if(linestream.fail()) throw IO::Error("Failed to parse species");
        Element &type = (*m.pte)[name];
        type.m = std::stod(mass);
        type.PWPP = pwpp;
    }
}

void parseCoordinates(std::string name, std::istream& file,
                      Molecule& m, IO::Parameter& p)
{
    auto &sys = std::get<NameList>(p.at("&SYSTEM"));
    auto dataentry = sys.find("nat");
    if (dataentry == sys.end()) throw IO::Error("nat not specified");
    size_t nat = std::stoul(dataentry->second);
    sys.erase(dataentry);
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
        line = IO::trim(line);
        while(line[0]=='!' || line[0]=='#'){
            std::getline(file, line);
            line = IO::trim(line);
        }
        std::stringstream linestream{line};
        linestream >> at.name;
        for(size_t i=0; i<3; ++i){
            linestream >> coord_expr;
            at.coord[i] = te_interp(coord_expr.c_str(), &err_pos);
            if(err_pos){
                throw IO::Error("Error parsing atom: "+coord_expr);
            }
        }
        if (linestream.fail()) {
            throw IO::Error{"Failed to parse atom"};
        }
        bool x{true},y{true},z{true};
        linestream >> x >> y >> z;
        at.properties->flags[AtomFlag::FixX] = !x;
        at.properties->flags[AtomFlag::FixY] = !y;
        at.properties->flags[AtomFlag::FixZ] = !z;
    }
}

void parseKPoints(std::string name, std::istream& file, Molecule& m)
{
    if (name.find("GAMMA") != name.npos) {
        return;
    } else if (name.find("AUTOMATIC") != name.npos) {
        std::string line;
        std::getline(file, line);
        line = IO::trim(line);
        while(line[0]=='!' || line[0]=='#'){
            std::getline(file, line);
            line = IO::trim(line);
        }
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

void parseCell(std::string name, std::istream& file, CellInp &cell)
{
    std::string line;
    for(size_t i=0; i<3; ++i){
        std::getline(file, line);
        line = IO::trim(line);
        while(line[0]=='!' || line[0]=='#'){
            std::getline(file, line);
            line = IO::trim(line);
        }
        std::stringstream linestream{line};
        linestream >> cell.cell[i][0] >> cell.cell[i][1] >> cell.cell[i][2];
        if (linestream.fail()) throw IO::Error("Failed to parse CELL_PARAMETERS");
    }
    if (name.find("BOHR") != name.npos) cell.fmt = CellInp::Bohr;
    else if (name.find("ANGSTROM") != name.npos) cell.fmt = CellInp::Angstrom;
    else cell.fmt = CellInp::Alat;
}

void createCell(Molecule &m, IO::Parameter &p, CellInp &cell)
{
    auto &s = m.getStep(0);
    auto &sys = std::get<NameList>(p.at("&SYSTEM"));
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
                    s.setCellDim(std::stod(celldm->second), CdmFmt::Bohr, scale);
                    sys.erase(celldm);
                } else if ((celldm == sys.end()) && (cellA != sys.end())) {
                    s.setCellDim(std::stod(cellA->second), CdmFmt::Angstrom, scale);
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
            axes[0] = stod(celldm->second);
            auto celldm2 = sys.find("celldm(2)");
            if(celldm2 != sys.end()) axes[1] = stod(celldm2->second);
            auto celldm3 = sys.find("celldm(3)");
            if(celldm3 != sys.end()) axes[2] = stod(celldm3->second);
            auto celldm4 = sys.find("celldm(4)");
            if(celldm4 != sys.end()) angles[0] = stod(celldm4->second);
            auto celldm5 = sys.find("celldm(5)");
            if(celldm5 != sys.end()) angles[1] = stod(celldm5->second);
            auto celldm6 = sys.find("celldm(6)");
            if(celldm6 != sys.end()) angles[2] = stod(celldm6->second);
        }else{
            axes[0] = stod(cellA->second);
            auto cellB = sys.find("B");
            if(cellB != sys.end()) axes[1] = stod(cellB->second)/axes[0];
            auto cellC = sys.find("C");
            if(cellC != sys.end()) axes[2] = stod(cellC->second)/axes[0];
            auto cosAB = sys.find("cosAB");
            if(cosAB != sys.end()) angles[0] = stod(cosAB->second);
            auto cosAC = sys.find("cosAC");
            if(cosAC != sys.end()) angles[1] = stod(cosAC->second);
            auto cosBC = sys.find("cosBC");
            if(cosBC != sys.end()) angles[2] = stod(cosBC->second);
        }
        s.setCellDim(axes[0], cdmFmt, s.getFmt() >= AtomFmt::Crystal);
        auto cell = makeBravais(ibrav, axes, angles);
        s.setCellVec(cell, s.getFmt() == AtomFmt::Crystal);
    }
}

void parseCard(std::string name, std::istream& file,
               Molecule& m, IO::Parameter& p,
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

IO::Data PWInpParser(const std::string& name, std::istream &file)
{
    IO::Data d{};
    Molecule &m = d.mol;
    m.setName(name);
    m.newStep();
    d.param = makeParam();
    auto &p = *d.param;
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

bool PWInpWriter(const Molecule& m, std::ostream &file,
                 const std::optional<IO::Parameter>& p,
                 const std::optional<IO::Preset>& c,
                 size_t index)
{
    if(!p || p->getFmt() != &IO::PWInput){
        throw IO::Error("PWI-Writer needs PWScf parameter set");
    }
    if(!c || c->getFmt() != &IO::PWInput){
        throw IO::Error("PWI-Writer needs suitable IO preset");
    }
    const auto &atfmt = std::get<NamedEnum>(c->at("atoms").first);
    const auto &cellfmt = std::get<NamedEnum>(c->at("atoms").first);
    const auto& s = (atfmt.name() == "Active") ?
        static_cast<const StepConst<Step::source>&>(m.getStep(index)) : // use active fmt
        m.getStep(index).asFmt(static_cast<AtomFmt>(atfmt.value())); // use explicit fmt
    const auto& PPPrefix = std::get<std::string>(p->at(""));
    const auto& PPSuffix = std::get<std::string>(p->at("PPSuffix"));
    const auto& control = std::get<NameList>(p->at("&CONTROL"));
    std::vector<std::string>
            outNL = {"&CONTROL", "&SYSTEM", "&ELECTRONS"};
    auto calc = control.find("calculation");
    if(calc != control.end()){
        if(calc->second == "'vc-relax'"){
            outNL.push_back("&IONS");
            outNL.push_back("&CELL");
        }else if(calc->second == "'relax'"){
            outNL.push_back("&IONS");
        }
    }
    for(auto &name: outNL){
        file << name << '\n';
        if(name == "&SYSTEM"){
            file << " ibrav = 0\n";
            file << " nat = " << s.getNat() << '\n';
            file << " ntyp = " << s.getNtyp() << '\n';
            auto cell_fmt = (cellfmt.name() == "Active") ?
                        ((s.getFmt() == AtomFmt::Angstrom) ?
                             CdmFmt::Angstrom : CdmFmt::Bohr) : // match coordinates
                        static_cast<CdmFmt>(cellfmt.value()); // use explicit
            if(cell_fmt == CdmFmt::Bohr){
                file << " celldm(1) = " << s.getCellDim(cell_fmt) << '\n';
            }else{
                file << " A = " << s.getCellDim(cell_fmt) << '\n';
            }
        }
        for(auto& e: std::get<NameList>(p->at(name))){
            file << ' ' << e.first << " = " << e.second << '\n';
        }
        file << "/\n\n";
    }
    file << "ATOMIC_SPECIES\n"
         << std::fixed << std::setprecision(5);
    for(auto &t: s.getTypes()){
        auto e = (*s.pte)[t];
        file << std::left << std::setw(3) << t << ' '
             << std::right << std::setw(9) << e.m << ' ';
        if(e.PWPP.empty()){
            file << PPPrefix << t << PPSuffix << '\n';
        }else{
            file << e.PWPP << '\n';
        }
    }
    const std::array<std::string, 4> fmt2str = {{"bohr", "angstrom", "crystal", "alat"}};
    file << "\nATOMIC_POSITIONS " << fmt2str[static_cast<size_t>(s.getFmt())] << '\n'
         << std::fixed << std::setprecision(5);
    AtomFlags fixComp{};
    fixComp[AtomFlag::FixX] = true;
    fixComp[AtomFlag::FixY] = true;
    fixComp[AtomFlag::FixZ] = true;
    for (const auto& at: s) {
        file << std::left << std::setw(3) << at.name
             << std::right << std::setprecision(8)
             << ' ' << std::setw(12) << at.coord[0]
             << ' ' << std::setw(12) << at.coord[1]
             << ' ' << std::setw(12) << at.coord[2];
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
    file << "\nCELL_PARAMETERS alat\n" << std::fixed << std::setprecision(8);
    for(auto &v: s.getCellVec()){
        file << std::setw(12) << v[0] << ' '
             << std::setw(12) << v[1] << ' '
             << std::setw(12) << v[2] << '\n';
    }
    return true;
}

const IO::Plugin IO::PWInput =
{
    "PWScf Input File",
    "pwi",
    "pwi",
    &PWInpParser,
    &PWInpWriter,
    &makeParam,
    &makePreset
};
