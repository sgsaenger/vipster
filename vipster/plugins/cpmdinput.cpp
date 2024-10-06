#include "cpmdinput.h"
#include "../ioutil.h"
#include <cctype>

using namespace Vipster;

using Section = std::vector<std::string>;
static const std::string sections[]{
    "&INFO", "&CPMD", "&SYSTEM", "&PIMD",
    "&PATH", "&PTDDFT", "&ATOMS", "&DFT",
    "&PROP", "&RESP", "&LINRES", "&TDDFT",
    "&HARDNESS", "&CLASSIC", "&EXTE", "&VDW",
    "&QMMM",
};

static const std::string SectionDoc{
    "Raw-text section that will be written to input file. "
    "See CPMD documentation for valid contents."
};

static Parameter makeParam()
{
    return {&Plugins::CPInput, {
            {"&INFO", {Section{}, SectionDoc}},
            {"&CPMD", {Section{}, SectionDoc}},
            {"&SYSTEM", {Section{}, SectionDoc+
                "\nVipster automatically fills in supported information:\n"
                "CELL VECTORS, KPOINTS and ANGSTROM/SCALE keywords according to I/O preset."}},
            {"&PIMD", {Section{}, SectionDoc}},
            {"&PATH", {Section{}, SectionDoc}},
            {"&PTDDFT", {Section{}, SectionDoc}},
            {"&ATOMS", {Section{}, SectionDoc+
                "\nVipster automatically fills in supported information:\n"
                "Atomic positions, ISOTOPE and pseudo-potential information and "
                "'FIX ATOMS' or 'FIX COORDINATES' constraints."}},
            {"&DFT", {Section{}, SectionDoc}},
            {"&PROP", {Section{}, SectionDoc}},
            {"&RESP", {Section{}, SectionDoc}},
            {"&LINRES", {Section{}, SectionDoc}},
            {"&TDDFT", {Section{}, SectionDoc}},
            {"&HARDNESS", {Section{}, SectionDoc}},
            {"&CLASSIC", {Section{}, SectionDoc}},
            {"&EXTE", {Section{}, SectionDoc}},
            {"&VDW", {Section{}, SectionDoc}},
            {"&QMMM", {Section{}, SectionDoc}},
            {"PPPrefix", {std::string{},
                    "The default prefix for pseudo-potentials.\n"
                    "Will be prepended to the type name if the periodic table entry for this type "
                    "does not contain an explicit pseudo-potential name."}},
            {"PPSuffix", {std::string{},
                    "The default suffix for pseudo-potentials.\n"
                    "Will be appended to the type name if the periodic table entry for this type "
                    "does not contain an explicit pseudo-potential name."}},
            {"PPNonlocality", {std::string{},
                    "The default non-locality for pseudo-potentials.\n"
                    "Will be printed as-is if the periodic table entry for the current type "
                    "does not contain an explicit non-locality string."}},
        }};
}

static Preset makePreset()
{
    return {&Plugins::CPInput,
        {{"fmt", {NamedEnum{4, {"Crystal", "Alat", "Angstrom", "Bohr", "Active"}},
                    "Active: Use the current Step's active atom format\n"
                    "Else: Enforce the selected atom format"}}}};
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

static Mat makeCell(int ibrav, double b, double c,
                    double alpha, double beta, double gamma)
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
                Vec{-0.5, 0.5*sqrt(3.), 0},
                Vec{   0, 0,             c}};
    case 5:
        // trigonal
    {
        double tx = sqrt((1 -   alpha)/2.);
        double ty = sqrt((1 -   alpha)/6.);
        double tz = sqrt((1 + 2*alpha)/3.);
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
        return {Vec{ 0.5, -0.5, c*0.5},
                Vec{ 0.5,  0.5, c*0.5},
                Vec{-0.5, -0.5, c*0.5}};
    case 0:
        [[fallthrough]];
    case 8:
        // simple orthorhombic
        return {Vec{1, 0, 0},
                Vec{0, b, 0},
                Vec{0, 0, c}};
    case 12:
        return {Vec{      1,                        0, 0},
                Vec{b*alpha, b*sqrt(1.-alpha*alpha), 0},
                Vec{      0,                        0, c}};
    case 14:
    {
        double singam = sqrt(1 - gamma*gamma);
        return {Vec{      1,        0, 0},
                Vec{b*gamma, b*singam, 0},
                Vec{ c*beta, c*(alpha-beta*gamma)/singam,
                     c*sqrt(1 + 2*alpha*beta*gamma - alpha*alpha - beta*beta - gamma*gamma)
                            /singam}};
    }
    default:
        break;
    }
    throw IOError("CPMD Input: Invalid ibrav");
}

IOTuple CPInpParser(const std::string& name, std::istream &file){
    Molecule m{name, 1};
    m.name = name;
    Step &s = m.getStep(0);
    Parameter p = makeParam();
    std::string curName{};
    Section *curSection{nullptr};


    std::string buf, line;
    int ibrav{0};
    double cellDim{1}, b{1}, c{1}, alpha{0}, beta{0}, gamma{0};
    AtomFmt cf{AtomFmt::Bohr};
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
                curName = "";
                curSection = nullptr;
            }else{
                auto pos = std::find_if(p.begin(), p.end(),
                                        [&line](const Parameter::value_type &pair){
                                            return pair.first == line;
                                        });
                if(pos != p.end()){
                    curName = line;
                    curSection = &std::get<Section>(pos->second.first);
                }else{
                    throw IOError("CPMD Input: Unknown section: "+line);
                }
            }
        }else{
            if(!curSection){
                throw IOError("CPMD Input: Data outside of section");
            }
            bool parsed{false};
            if(curName == "&SYSTEM"){
                if(line.find("ANGSTROM") != line.npos) {
                    parsed = true;
                    cf = AtomFmt::Angstrom;
                    af = AtomFmt::Angstrom;
                }else if(line.find("KPOINTS") != line.npos){
                    parsed = true;
                    auto& kp = m.kpoints;
                    if(line.find("MONKHORST-PACK") != line.npos){
                        kp.active = KPoints::Fmt::MPG;
                        std::getline(file, buf);
                        std::stringstream{buf} >> kp.mpg.x >> kp.mpg.y >> kp.mpg.z;
                        size_t pos;
                        if((pos = line.find("SHIFT=")) != line.npos){
                            std::stringstream{line.substr(pos+6)}
                                >> kp.mpg.sx >> kp.mpg.sy >> kp.mpg.sz;
                        }
                    }else{
                        kp.active = KPoints::Fmt::Discrete;
                        if(line.find("SCALED") != line.npos){
                            kp.discrete.properties |=
                                    KPoints::Discrete::Properties::crystal;
                        }
                        if(line.find("BANDS") != line.npos){
                            kp.discrete.properties |=
                                    KPoints::Discrete::Properties::band;
                            auto isTerm = [](const KPoints::Discrete::Point& p1,
                                             const KPoints::Discrete::Point& p2){
                                return float_comp(p1.pos[0], 0.) && float_comp(p1.pos[1], 0.)
                                        && float_comp(p1.pos[2], 0.) && float_comp(p1.weight, 0.)
                                        && float_comp(p2.pos[0], 0.) && float_comp(p2.pos[1], 0.)
                                        && float_comp(p2.pos[2], 0.);
                            };
                            bool cont{true};
                            do{
                                std::getline(file, buf);
                                KPoints::Discrete::Point p1, p2;
                                std::stringstream{buf}
                                    >> p1.weight
                                    >> p1.pos[0] >> p1.pos[1] >> p1.pos[2]
                                    >> p2.pos[0] >> p2.pos[1] >> p2.pos[2];
                                if(isTerm(p1, p2)){
                                    cont = false;
                                }else{
                                    kp.discrete.kpoints.push_back(p1);
                                    kp.discrete.kpoints.push_back(p2);
                                }
                            }while(cont);
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
                        double fac;
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
            }else if(curName == "&ATOMS"){
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
                            name2 = name + std::to_string(++i);
                        }
                        name = name2;
                    }
                    types.push_back(name);
                    s.getPTE()[name].CPPP = CPPP;
                    std::getline(file, buf);
                    s.getPTE()[name].CPNL = trim(buf);
                    std::getline(file, buf);
                    size_t oldNat = s.getNat();
                    size_t nat = std::stoul(buf);
                    s.newAtoms(nat);
                    for(auto at = s.begin()+oldNat; at!=s.end(); ++at){
                        std::getline(file, buf);
                        std::stringstream linestream{buf};
                        at->name = name;
                        linestream >> at->coord[0] >> at->coord[1] >> at->coord[2];
                        if (linestream.fail()) {
                            throw IOError{"CPMD Input: Failed to parse atom"};
                        }
                    }
                }else if(line.find("CONSTRAINTS") != line.npos){
                    std::getline(file, buf);
                    parsed = true;
                    std::vector<std::string> other;
                    while(buf.find("END CONSTRAINTS") == buf.npos){
                        AtomProperties::Flags fixComp{};
                        fixComp[AtomProperties::FixX] = true;
                        fixComp[AtomProperties::FixY] = true;
                        fixComp[AtomProperties::FixZ] = true;
                        if(buf.find("FIX") != buf.npos){
                            if((buf.find("ALL") != buf.npos) ||
                               (buf.find("QM") != buf.npos)){
                                for(auto& at: s){
                                    at.properties->flags |= fixComp;
                                }
                            }else if(buf.find("ELEM") != buf.npos){
                                bool seq = buf.find("SEQ") != buf.npos;
                                std::getline(file, buf);
                                size_t Z;
                                auto it=s.begin();
                                auto it_end=s.end();
                                if(seq){
                                    size_t beg, end;
                                    std::stringstream{buf} >> Z >> beg >> end;
                                    it += beg-1;
                                    it_end = s.begin() + end;
                                }else{
                                    std::stringstream{buf} >> Z;
                                }
                                for(;it != it_end; ++it){
                                    if(it->type->Z == Z){
                                        it->properties->flags |= fixComp;
                                    }
                                }
                            }else if(buf.find("PPTY") != buf.npos){
                                bool seq = buf.find("SEQ") != buf.npos;
                                std::getline(file, buf);
                                size_t pp;
                                auto it=s.begin();
                                auto it_end=s.end();
                                if(seq){
                                    size_t beg, end;
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
                                size_t nat, idx;
                                file >> nat;
                                for(size_t i=0; i<nat; ++i){
                                    file >> idx;
                                    s[idx-1].properties->flags |= fixComp;
                                }
                            }else if(buf.find("COOR") != buf.npos){
                                size_t nat, idx;
                                file >> nat;
                                for(size_t i=0; i < nat; ++i){
                                    bool x{}, y{}, z{};
                                    file >> idx >> x >> y >> z;
                                    s[idx-1].properties->flags[AtomProperties::FixX] = !x;
                                    s[idx-1].properties->flags[AtomProperties::FixY] = !y;
                                    s[idx-1].properties->flags[AtomProperties::FixZ] = !z;
                                }
                            }else{
                                other.push_back(buf);
                            }
                        }else if(!buf.empty()){
                            other.push_back(buf);
                        }
                        std::getline(file, buf);
                    }
                    if(!other.empty()){
                        curSection->push_back("CONSTRAINTS");
                        curSection->insert(curSection->end(),
                                           other.begin(),
                                           other.end());
                        curSection->push_back("END CONSTRAINTS");
                    }
                }else if(line.find("ISOTOPE") != line.npos){
                    parsed = true;
                    for(auto& t: types){
                        std::getline(file, buf);
                        s.getPTE()[t].m = std::stod(buf);
                    }
                }
            }
            if(!parsed){
                curSection->push_back(buf);
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
    s.setFmt(af, false);
    if(scale){
        for(auto& at: s){
            at.coord[0] /= scaleVec[0];
            at.coord[1] /= scaleVec[1];
            at.coord[2] /= scaleVec[2];
        }
    }

    return {std::move(m), std::move(p), DataList{}};
}

bool CPInpWriter(const Molecule& m, std::ostream &file,
                 const std::optional<Parameter>& p,
                 const std::optional<Preset>& c,
                 size_t index)
{
    if(!p || p->getFmt() != &Plugins::CPInput){
        throw IOError("CPMD Input: writer needs CPMD parameter set");
    }
    if(!c || c->getFmt() != &Plugins::CPInput){
        throw IOError("CPMD Input: writer needs suitable IO preset");
    }
    const auto &atfmt = std::get<NamedEnum>(c->at("fmt").first);
    const auto& s =  m.getStep(index).asFmt((atfmt.name() == "Active") ?
                                            m.getStep(index).getFmt() : // use active format
                                            static_cast<AtomFmt>(atfmt.value()-2)); // use explicit format
    auto cf = (s.getFmt() == AtomFmt::Angstrom) ? AtomFmt::Angstrom : AtomFmt::Bohr;
    const auto& PPPrefix = std::get<std::string>(p->at("PPPrefix").first);
    const auto& PPSuffix = std::get<std::string>(p->at("PPSuffix").first);
    const auto& PPNonlocality = std::get<std::string>(p->at("PPNonlocality").first);
    for(const auto& name: sections){
        file << name << '\n';
        const auto& curSection = std::get<Section>(p->at(name).first);
        if(name == "&ATOMS"){
            std::map<std::string, std::vector<size_t>> types;
            for(auto it=s.begin(); it!=s.end(); ++it){
                types[it->name].push_back(it.getIdx());
            }
            AtomProperties::Flags fixComp{};
            fixComp[AtomProperties::FixX] = true;
            fixComp[AtomProperties::FixY] = true;
            fixComp[AtomProperties::FixZ] = true;
            size_t count{0};
            std::vector<size_t> fixAtom;
            std::vector<std::pair<size_t, AtomProperties::Flags>> fixCoord;
            std::vector<double> masses;
            for(const auto& pair: types){
                const auto& pE = m.getPTE().at(pair.first);
                masses.push_back(pE.m);
                if(!pE.CPPP.empty()){
                    file << '*' << pE.CPPP << '\n';
                }else{
                    file << '*' << trim(PPPrefix) << trim(pair.first)
                         << trim(PPSuffix) << '\n';
                }
                if(!pE.CPNL.empty()){
                    file << "  " << pE.CPNL << '\n';
                }else{
                    file << "  " << PPNonlocality << '\n';
                }
                file << "  " << pair.second.size() << '\n'
                     << std::fixed << std::setprecision(5);
                auto it = s.begin();
                for(const auto& i: pair.second){
                    it += i-it.getIdx();
                    count += 1;
                    file << ' '
                         << ' ' << std::right << std::setw(10) << it->coord[0]
                         << ' ' << std::right << std::setw(10) << it->coord[1]
                         << ' ' << std::right << std::setw(10) << it->coord[2] << '\n';
                    const auto fix = it->properties->flags & fixComp;
                    if(fix == fixComp){
                        fixAtom.push_back(count);
                    }else if(fix.any()){
                        fixCoord.emplace_back(count, fix);
                    }
                }
            }
            if(fixAtom.empty() && fixCoord.empty()){
                for(const auto& line: curSection){
                    file << line << '\n';
                }
            }else{
                size_t otherConstraints{curSection.size()};
                for(size_t i=0; i<curSection.size(); ++i){
                    const auto& line = curSection[i];
                    file << line << '\n';
                    if(line.find("CONSTRAINTS") != line.npos){
                        otherConstraints = i;
                        break;
                    }
                }
                if(otherConstraints == curSection.size()){
                    file << "  CONSTRAINTS\n";
                }
                if(!fixAtom.empty()){
                    file << "  FIX ATOMS\n  "
                         << fixAtom.size() << ' ';
                    for(const auto& idx: fixAtom){
                        file << idx << ' ';
                    }
                    file << "\n\n";
                }
                if(!fixCoord.empty()){
                    file << "  FIX COORDINATES\n  "
                         << fixCoord.size() << "\n  ";
                    for(const auto& pair: fixCoord){
                        file << pair.first << ' '
                             << pair.second[AtomProperties::FixX] << ' '
                             << pair.second[AtomProperties::FixY] << ' '
                             << pair.second[AtomProperties::FixZ] << "\n  ";
                    }
                    file << '\n';
                }
                if(otherConstraints != curSection.size()){
                    for(size_t i=otherConstraints+1; i<curSection.size(); ++i){
                        file << curSection[i] << '\n';
                    }
                }else{
                    file << "  END CONSTRAINTS\n";
                }
            }
            file << "  ISOTOPE\n";
            for(const auto& m: masses){
                file << "  " << m << '\n';
            }
            file << "&END\n\n";
        }else if(name == "&SYSTEM"){
            switch(s.getFmt()){
            case AtomFmt::Angstrom:
                file << "  ANGSTROM\n";
                break;
            case AtomFmt::Alat:
                file << "  SCALE CARTESIAN\n";
                break;
            case AtomFmt::Crystal:
                file << "  SCALE\n";
                break;
            default:
                break;
            }
            Mat tmpvec = s.getCellVec() * s.getCellDim(cf);
            file << "  CELL VECTORS\n" << std::fixed << std::setprecision(5)
                 << "  " << tmpvec[0][0] << ' ' << tmpvec[0][1] << ' ' << tmpvec[0][2] << '\n'
                 << "  " << tmpvec[1][0] << ' ' << tmpvec[1][1] << ' ' << tmpvec[1][2] << '\n'
                 << "  " << tmpvec[2][0] << ' ' << tmpvec[2][1] << ' ' << tmpvec[2][2] << '\n';
            const auto kp = m.kpoints;
            if(kp.active == KPoints::Fmt::MPG){
                const auto& mpg = kp.mpg;
                file << "  KPOINTS MONKHORST-PACK";
                if(float_comp(mpg.sx, 0.) || float_comp(mpg.sy, 0.) || float_comp(mpg.sz, 0.)){
                    file << " SHIFT=" << mpg.sx << ' ' << mpg.sy << ' ' << mpg.sz;
                }
                file << "\n  " << mpg.x << ' ' << mpg.y << ' ' << mpg.z << '\n';
            }else if(kp.active == KPoints::Fmt::Discrete){
                const auto& disc = kp.discrete;
                file << "  KPOINTS" << std::defaultfloat;
                if(disc.properties & KPoints::Discrete::crystal){
                    file << " SCALED";
                }
                if(disc.properties & KPoints::Discrete::band){
                    file << " BANDS\n";
                    if(disc.kpoints.size()%2){
                        throw IOError("CPMD Input: for BANDS, number of K-Points needs to be even");
                    }
                    for(auto it=disc.kpoints.begin(); it<disc.kpoints.end(); it+=2){
                        const auto& k1 = *it;
                        const auto& k2 = *(it+1);
                        file << "  " << std::floor(k1.weight)
                             << ' ' << k1.pos[0] << ' ' << k1.pos[1] << ' ' << k1.pos[2]
                             << ' ' << k2.pos[0] << ' ' << k2.pos[1] << ' ' << k2.pos[2];
                    }
                    file << "0 0. 0. 0. 0. 0. 0.\n";
                }else{
                    file << "\n  " << disc.kpoints.size() << '\n';
                    for(const auto& k: disc.kpoints){
                        file << "  " << k.pos[0] << ' ' << k.pos[1]
                             << ' ' << k.pos[2] << ' ' << k.weight << '\n';
                    }
                }
            }
            for(const auto& line: curSection){
                file << line << '\n';
            }
            file << "&END\n\n";
        }else if(name == "&CPMD" ||
                 !curSection.empty()){
            for(const auto& line: curSection){
                file << line << '\n';
            }
            file << "&END\n\n";
        }
    }
    return true;
}

const Plugin Plugins::CPInput =
{
    "CPMD Input File",
    "cpi",
    "cpi",
    &CPInpParser,
    &CPInpWriter,
    &makeParam,
    &makePreset
};
