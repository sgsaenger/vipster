#include "plugin.h"

using namespace Vipster;

static std::string trim(const std::string& str)
{
    const std::string whitespace{" \t"};
    const auto strBegin = str.find_first_not_of(whitespace);
    if(strBegin == std::string::npos){
        return "";
    }
    const auto strEnd = str.find_last_not_of(whitespace);
    return str.substr(strBegin, strEnd-strBegin+1);
}

const std::map<std::string, int> str2ibrav{
    {"ISOLATED", 0},
    {"CUBIC", 1},
    {"FACE CENTERED CUBIC", 2},
    {"FCC", 2},
    {"BODY CENTERED CUBIC", 3},
    {"BCC", 3},
    {"HEXAGONAL", 4},
    {"TRIGONAL", 5},
    {"RHOMBOHEDRAL", 5},
    {"TETRAGONAL", 6},
    {"BODY CENTERED TETRAGONAL", 7},
    {"BCT", 7},
    {"ORTHORHOMBIC", 8},
    {"MONOCLINIC", 12},
    {"TRICLINIC", 14}
};

const std::map<std::string, IO::CPParam::Section IO::CPParam::*> IO::CPParam::str2section{
    {"&INFO", &IO::CPParam::info},
    {"&CPMD", &IO::CPParam::cpmd},
    {"&SYSTEM", &IO::CPParam::system},
    {"&PIMD", &IO::CPParam::pimd},
    {"&PATH", &IO::CPParam::path},
    {"&PTDDFT", &IO::CPParam::ptddft},
    {"&ATOMS", &IO::CPParam::atoms},
    {"&DFT", &IO::CPParam::dft},
    {"&PROP", &IO::CPParam::prop},
    {"&RESP", &IO::CPParam::resp},
    {"&LINRES", &IO::CPParam::linres},
    {"&TDDFT", &IO::CPParam::tddft},
    {"&HARDNESS", &IO::CPParam::hardness},
    {"&CLASSIC", &IO::CPParam::classic},
    {"&EXTE", &IO::CPParam::exte},
    {"&VDW", &IO::CPParam::vdw},
    {"&QMMM", &IO::CPParam::qmmm},
};

static Mat makeCell(int ibrav, float b, float c,
                    float alpha, float beta, float gamma)
{
    switch(ibrav){
    case 1:
        // simple cubic
        return {Vec{1, 0, 0},
                Vec{0, 1, 0},
                Vec{0, 0, 1}};
    case 2:
        // fcc
        return {Vec{-0.5,    0,  0.5},
                Vec{   0,  0.5,  0.5},
                Vec{-0.5,  0.5,    0}};
    case 3:
        // bcc
        return {Vec{ 0.5,  0.5,  0.5},
                Vec{-0.5,  0.5,  0.5},
                Vec{-0.5, -0.5,  0.5}};
    case 4:
        // hexagonal
        return {Vec{   1, 0.0,           0},
                Vec{-0.5, 0.5f*sqrtf(3), 0},
                Vec{   0, 0,             c}};
    case 5:
        // trigonal
    {
        float tx = sqrtf((1 -   alpha)/2.f);
        float ty = sqrtf((1 -   alpha)/6.f);
        float tz = sqrtf((1 + 2*alpha)/3.f);
        return {Vec{ tx,  -ty, tz},
                Vec{  0, 2*ty, tz},
                Vec{-tx,  -ty, tz}};
    }
    case 6:
        // simple tetragonal
        return {Vec{1, 0, 0},
                Vec{0, 1, 0},
                Vec{0, 0, c}};
    case 7:
        // body centered tetragonal
        return {Vec{ 0.5, -0.5, c*0.5f},
                Vec{ 0.5,  0.5, c*0.5f},
                Vec{-0.5, -0.5, c*0.5f}};
    case 0:
        [[fallthrough]];
    case 8:
        // simple orthorhombic
        return {Vec{1, 0, 0},
                Vec{0, b, 0},
                Vec{0, 0, c}};
    case 12:
        return {Vec{      1,                        0, 0},
                Vec{b*alpha, b*sqrtf(1.f-alpha*alpha), 0},
                Vec{      0,                        0, c}};
    case 14:
    {
        float singam = sqrtf(1 - gamma*gamma);
        return {Vec{      1,        0, 0},
                Vec{b*gamma, b*singam, 0},
                Vec{ c*beta, c*(alpha-beta*gamma)/singam,
                     c*sqrtf(1 + 2*alpha*beta*gamma - alpha*alpha - beta*beta - gamma*gamma)
                            /singam}};
    }
    default:
        break;
    }
    throw IO::Error("Invalid ibrav");
}

IO::Data CPInpParser(const std::string& name, std::ifstream &file){
    IO::Data d{};
    d.fmt = IOFmt::CPI;
    Molecule &m = d.mol;
    m.setName(name);
    StepProper &s = m.newStep();
    d.param = std::make_unique<IO::CPParam>(name);
    IO::CPParam &p = *static_cast<IO::CPParam*>(d.param.get());
    IO::CPParam::Section IO::CPParam::* curSection{nullptr};

    std::string buf, line;
    int ibrav{0};
    float cellDim{1}, b{1}, c{1}, alpha{0}, beta{0}, gamma{0};
    CdmFmt cf{CdmFmt::Bohr};
    AtomFmt af{AtomFmt::Bohr};
    std::vector<std::string> types;
    bool scale{false};
    Vec scaleVec{1,1,1};
    Mat cellVec{Vec{1,0,0},Vec{0,1,0},Vec{0,0,1}};
    while (std::getline(file, buf)) {
        line = trim(buf);
        if(line.empty() || line[0] == '!') continue;
        for (auto &c: line) c = static_cast<char>(std::toupper(c));
        if(line[0] == '&'){
            if(line.find("&END") != line.npos){
                curSection = nullptr;
            }else{
                curSection = IO::CPParam::str2section.at(line);
            }
        }else{
            if(!curSection){
                throw IO::Error("Data outside of section");
            }
            bool parsed{false};
            if(curSection == &IO::CPParam::system){
                if(line.find("ANGSTROM") != line.npos) {
                    parsed = true;
                    cf = CdmFmt::Angstrom;
                    af = AtomFmt::Angstrom;
                }else if(line.find("KPOINTS") != line.npos){
                    parsed = true;
                    auto& kp = m.getKPoints();
                    if(line.find("MONKHORST-PACK") != line.npos){
                        kp.active = KPointFmt::MPG;
                        std::getline(file, buf);
                        std::stringstream{buf} >> kp.mpg.x >> kp.mpg.y >> kp.mpg.z;
                        size_t pos;
                        if((pos = line.find("SHIFT=")) != line.npos){
                            std::stringstream{line.substr(pos+6)}
                                >> kp.mpg.sx >> kp.mpg.sy >> kp.mpg.sz;
                        }
                    }else{
                        kp.active = KPointFmt::Discrete;
                        if(line.find("SCALED") != line.npos){
                            kp.discrete.properties |=
                                    KPoints::Discrete::Properties::crystal;
                        }
                        if(line.find("BANDS") != line.npos){
                            kp.discrete.properties |=
                                    KPoints::Discrete::Properties::band;
                            auto isTerm = [](const DiscreteKPoint& p){
                                return float_comp(p.pos[0], 0) && float_comp(p.pos[1], 0)
                                        && float_comp(p.pos[2], 0) && float_comp(p.weight, 0);
                            };
                            do{
                                std::getline(file, buf);
                                DiscreteKPoint p1, p2;
                                std::stringstream{buf}
                                    >> p1.weight
                                    >> p1.pos[0] >> p1.pos[1] >> p1.pos[2]
                                    >> p2.pos[0] >> p2.pos[1] >> p2.pos[2];
                            }while(!isTerm(kp.discrete.kpoints.back()));
                        }else{
                            std::getline(file, buf);
                            size_t nk = std::stoul(buf);
                            kp.discrete.kpoints.resize(nk);
                            for(auto& p: kp.discrete.kpoints){
                                std::getline(file, buf);
                                std::stringstream{buf}
                                    >> p.pos[0] >> p.pos[1] >> p.pos[2] >> p.weight;
                            }
                        }
                    }
                }else if(line.find("SCALE") != line.npos){
                    parsed = true;
                    if(line.find("CARTESIAN") != line.npos){
                        af = AtomFmt::Alat;
                    }else{
                        af = AtomFmt::Crystal;
                    }
                    size_t pos;
                    if((pos = line.find("S=")) != line.npos){
                        scale = true;
                        float fac;
                        std::stringstream{line.substr(pos+2)} >> fac;
                        scaleVec = Vec{fac,fac,fac};
                    }
                    if((pos = line.find("SX=")) != line.npos){
                        scale = true;
                        std::stringstream{line.substr(pos+3)} >> scaleVec[0];
                    }
                    if((pos = line.find("SY=")) != line.npos){
                        scale = true;
                        std::stringstream{line.substr(pos+3)} >> scaleVec[1];
                    }
                    if((pos = line.find("SZ=")) != line.npos){
                        scale = true;
                        std::stringstream{line.substr(pos+3)} >> scaleVec[2];
                    }
                }else if(line.find("SYMMETRY") != line.npos){
                    parsed = true;
                    std::getline(file, buf);
                    std::stringstream ss{buf};
                    int tmp;
                    ss >> tmp;
                    if(ss.fail()){
                        buf = trim(buf);
                        for (auto &c: buf) c = static_cast<char>(std::toupper(c));
                        ibrav = str2ibrav.at(buf);
                    }else{
                        ibrav = tmp;
                    }
                }else if(line.find("CELL") != line.npos){
                    parsed = true;
                    if(line.find("VECTORS") != line.npos){
                        ibrav = -1;
                        file >> cellVec[0][0] >> cellVec[0][1] >> cellVec[0][2]
                             >> cellVec[1][0] >> cellVec[1][1] >> cellVec[1][2]
                             >> cellVec[2][0] >> cellVec[2][1] >> cellVec[2][2];
                    }else{
                        std::getline(file, buf);
                        std::stringstream{buf} >> cellDim >> b >> c >> alpha >> beta >> gamma;
                        if(line.find("ABSOLUTE") != line.npos){
                            b /= cellDim;
                            c /= cellDim;
                        }
                        if(line.find("DEGREE") != line.npos){
                            alpha = std::cos(deg2rad * alpha);
                            beta = std::cos(deg2rad * beta);
                            gamma = std::cos(deg2rad * gamma);
                        }
                    }
                }
            }else if(curSection == &IO::CPParam::atoms){
                if(line[0] == '*'){
                    parsed = true;
                    line = trim(buf); // reset to regain capitalization
                    std::string CPPP = line.substr(1);
                    auto pos = CPPP.find_first_of(". ");
                    std::string name = CPPP.substr(0, pos);
                    auto tmp = s.getTypes();
                    if(tmp.find(name) != tmp.end()){
                        int i = 1;
                        std::string name2 = name + std::to_string(i);
                        while(tmp.find(name2) != tmp.end()){
                            std::string name2 = name + std::to_string(++i);
                        }
                        name = name2;
                    }
                    types.push_back(name);
                    (*s.pse)[name].CPPP = CPPP;
                    std::getline(file, buf);
                    (*s.pse)[name].CPNL = trim(buf);
                    std::getline(file, buf);
                    size_t nat = std::stoul(buf);
                    s.newAtoms(nat);
                    for(auto& at: s){
                        std::getline(file, buf);
                        std::stringstream linestream{buf};
                        at.name = name;
                        linestream >> at.coord[0] >> at.coord[1] >> at.coord[2];
                        if (linestream.fail()) {
                            throw IO::Error{"Failed to parse atom"};
                        }
                    }
                }else if(line.find("CONSTRAINTS") != line.npos){
                    parsed = true;
                    std::getline(file, buf);
                    std::vector<std::string> other;
                    while(buf.find("END CONSTRAINTS") == buf.npos){
                        std::getline(file, buf);
                        AtomFlags fixComp{};
                        fixComp[AtomFlag::FixX] = true;
                        fixComp[AtomFlag::FixY] = true;
                        fixComp[AtomFlag::FixZ] = true;
                        if(buf.find("FIX") != buf.npos){
                            if((buf.find("ALL") != buf.npos) ||
                               (buf.find("QM") != buf.npos)){
                                for(auto& at: s){
                                    at.properties->flags |= fixComp;
                                }
                            }else if(buf.find("ELEM") != buf.npos){
                                bool seq = buf.find("SEQ") != buf.npos;
                                std::getline(file, buf);
                                size_t Z, beg, end;
                                auto it=s.begin();
                                auto it_end=s.end();
                                if(seq){
                                    std::stringstream{buf} >> Z >> beg >> end;
                                    it += beg-1;
                                    it_end = s.begin() + end;
                                }else{
                                    std::stringstream{buf} >> Z;
                                }
                                for(;it != it_end; ++it){
                                    if(it->pse->Z == Z){
                                        it->properties->flags |= fixComp;
                                    }
                                }
                            }else if(buf.find("PPTY") != buf.npos){
                                bool seq = buf.find("SEQ") != buf.npos;
                                std::getline(file, buf);
                                size_t pp, beg, end;
                                auto it=s.begin();
                                auto it_end=s.end();
                                if(seq){
                                    std::stringstream{buf} >> pp >> beg >> end;
                                    it += beg-1;
                                    it_end = s.begin() + end;
                                }else{
                                    std::stringstream{buf} >> pp;
                                }
                                const auto& type = types[pp];
                                for(;it != it_end; ++it){
                                    if(it->name == type){
                                        it->properties->flags |= fixComp;
                                    }
                                }
                            }else if(buf.find("SEQ") != buf.npos){
                                std::getline(file, buf);
                                size_t beg, end;
                                std::stringstream{buf} >> beg >> end;
                                for(auto it=s.begin()+(beg-1); it!=s.begin()+end; ++it){
                                    it->properties->flags |= fixComp;
                                }
                            }else if(buf.find("ATOM") != buf.npos){
                                std::getline(file, buf);
                                std::stringstream ss{buf};
                                size_t nat, idx;
                                ss >> nat;
                                for(size_t i=0; i<nat; ++i){
                                    ss >> idx;
                                    s[idx].properties->flags |= fixComp;
                                }
                            }else if(buf.find("COOR") != buf.npos){
                                std::getline(file, buf);
                                std::stringstream ss{buf};
                                size_t nat, idx;
                                ss >> nat;
                                for(size_t i=0; i < nat; ++i){
                                    bool x{}, y{}, z{};
                                    ss >> idx >> x >> y >> z;
                                    s[idx].properties->flags[AtomFlag::FixX] = x;
                                    s[idx].properties->flags[AtomFlag::FixY] = y;
                                    s[idx].properties->flags[AtomFlag::FixZ] = z;
                                }
                            }else{
                                other.push_back(buf);
                            }
                        }else{
                            other.push_back(buf);
                        }
                        if(!buf.empty()){
                            (p.*curSection).push_back("CONSTRAINTS");
                            (p.*curSection).insert((p.*curSection).end(),
                                                   other.begin(),
                                                   other.end());
                            (p.*curSection).push_back("END CONSTRAINTS");
                        }
                    }
                }else if(line.find("ISOTOPE") != line.npos){
                    parsed = true;
                    for(auto& t: types){
                        std::getline(file, buf);
                        (*s.pse)[t].m = std::stof(buf);
                    }
                }
            }
            if(!parsed){
                (p.*curSection).push_back(buf);
            }
        }
    }
    // create cellVec if needed, apply
    if(ibrav != -1){
        cellVec = makeCell(ibrav, b, c, alpha, beta, gamma);
    }
    s.setCellVec(cellVec);

    // apply scaling
    s.setCellDim(cellDim, cf);
    s.setFmt(af, true);
    if(scale){
        for(auto& at: s){
            at.coord[0] /= scaleVec[0];
            at.coord[1] /= scaleVec[1];
            at.coord[2] /= scaleVec[2];
        }
    }

    return d;
}

bool CPInpWriter(const Molecule& m, std::ofstream &file,
                 const BaseParam *const p,
                 const BaseConfig *const,
                 IO::State state)
{
    const auto& s = m.getStep(state.index);
    const auto *pp = dynamic_cast<const IO::CPParam*>(p);
    if(!pp) throw IO::Error("CPI-Writer needs CPMD parameter set");
    return true;
}

const IO::Plugin IO::CPInput =
{
    "CPMD Input File",
    "cpi",
    "cpi",
    IO::Plugin::Param,
    &CPInpParser,
    &CPInpWriter
};
