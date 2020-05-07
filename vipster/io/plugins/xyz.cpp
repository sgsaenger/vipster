#include "xyz.h"

#include <sstream>
#include <iomanip>

using namespace Vipster;

static IO::Preset makePreset()
{
    return {&IO::XYZ,
        {{"filemode", {NamedEnum{0, {"Step", "Trajec", "Cell"}},
                    "Step: Write currently active Step\n"
                    "Trajec: Write complete trajectory of active Molecule\n"
                    "Cell: Write currently active Step and append the unit cell"}},
         {"atomdata", {NamedEnum{0, {"None", "Charge", "Forces"}},
                    "None: Write only coordinates (standard xyz)\n"
                    "Charge: Append a column containing the atomic charges\n"
                    "Forces: Append three columns containing the atomic forces"}}}};
}


IO::Data XYZParser(const std::string& name, std::istream &file)
{
    IO::Data data{};
    Molecule &m = data.mol;
    m.name = name;

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream natline{line};
        size_t nat;
        natline >> nat;
        if (natline.fail()) {
            if(m.getNstep() == 1u){
                double f1{};
                std::getline(file, line);
                natline = std::stringstream{line};
                if(!(natline >> f1).fail()){
                    double f2{}, f3{};
                    if(!(natline >> f2 >> f3).fail()){
                        // this could be a cell-vector,
                        // expect next two lines to follow suit
                        auto tmp = Mat{Vec{f1,f2,f3},Vec{},Vec{}};
                        std::getline(file, line);
                        natline = std::stringstream{line};
                        natline >> tmp[1][0] >> tmp[1][1] >> tmp[1][2];
                        if(natline.fail()){
                            throw IO::Error("XYZ: Non-standard data after XYZ-file");
                        }
                        std::getline(file, line);
                        natline = std::stringstream{line};
                        natline >> tmp[2][0] >> tmp[2][1] >> tmp[2][2];
                        if(natline.fail()){
                            throw IO::Error("XYZ: Non-standard data after XYZ-file");
                        }
                        m.getStep(0).setCellVec(tmp);
                        m.getStep(0).setCellDim(1, AtomFmt::Angstrom);
                        break;
                    }else if((floor(f1) == f1) && (f1>0)){
                        // found nat, continue parsing trajectory
                        nat = static_cast<size_t>(f1);
                    }else{
                        throw IO::Error("XYZ: Failed to parse nat");
                    }
                }
            }else{
                while ((natline>>nat).fail()) {
                    std::getline(file, line);
                    if (file.eof()) {
                        if(m.getNstep() == 0u){
                            throw IO::Error("XYZ: Failed to parse nat");
                        }else{
                            throw IO::Error("XYZ: Non-standard data after XYZ-file");
                        }
                    }
                    natline = std::stringstream{line};
                }
            }
        }
        Step &sp = m.newStep();
        sp.setFmt(AtomFmt::Angstrom);
        sp.enableCell(false);
        sp.newAtoms(nat);
        std::getline(file, line);
        sp.setComment(line);
        for (auto& at: sp) {
            std::getline(file, line);
            std::stringstream atline{line};
            double f1{};
            atline >> at.name >> at.coord[0] >> at.coord[1] >> at.coord[2];
            if (atline.fail()) {
                throw IO::Error("XYZ: failed to parse atom");
            }
            if(!(atline >> f1).fail()){
                double f2{}, f3{};
                if(!(atline >> f2 >> f3).fail()){
                    at.properties->forces = {f1, f2, f3};
                }else{
                    at.properties->charge = f1;
                }
            }
        }
    }
    return data;
}

bool XYZWriter(const Molecule& m, std::ostream &file,
               const std::optional<IO::Parameter>&,
               const std::optional<IO::Preset>& c,
               size_t index)
{
    if(!c || c->getFmt() != &IO::XYZ){
        throw IO::Error("XYZ: writer needs suitable IO preset");
    }
    auto cc = *c;
    const auto& s = m.getStep(index).asFmt(AtomFmt::Angstrom);
    auto stepWriter = [&file, cc](const auto& s){
        file << s.getNat() << '\n';
        file << s.getComment() << '\n';
        file << std::fixed << std::setprecision(5);
        switch(std::get<NamedEnum>(cc.at("atomdata").first)){
        case 0: // None
            for(const auto& at: s){
                file << std::left << std::setw(3) << at.name << " "
                     << std::right << std::setw(10) << at.coord[0] << " "
                     << std::right << std::setw(10) << at.coord[1] << " "
                     << std::right << std::setw(10) << at.coord[2] << '\n';
            }
            break;
        case 1: // charge
            for(const auto& at: s){
                file << std::left << std::setw(3) << at.name << " "
                     << std::right << std::setw(10) << at.coord[0] << " "
                     << std::right << std::setw(10) << at.coord[1] << " "
                     << std::right << std::setw(10) << at.coord[2] << " "
                     << std::right << std::setw(10) << at.properties->charge << '\n';
            }
            break;
        case 2: // forces
            for(const auto& at: s){
                file << std::left << std::setw(3) << at.name << " "
                     << std::right << std::setw(10) << at.coord[0] << " "
                     << std::right << std::setw(10) << at.coord[1] << " "
                     << std::right << std::setw(10) << at.coord[2] << " "
                     << std::right << std::setw(10) << at.properties->forces[0] << " "
                     << std::right << std::setw(10) << at.properties->forces[1] << " "
                     << std::right << std::setw(10) << at.properties->forces[2] << '\n';
            }
            break;
        default:
            throw IO::Error{"XYZ: invalid atomdata setting"};
        }
    };
    switch(std::get<NamedEnum>(cc.at("filemode").first)){
    case 0: // Step
        stepWriter(s);
        break;
    case 1: // Trajec
        for(const auto& step: m.getSteps()){
            stepWriter(step.asFmt(AtomFmt::Angstrom));
        }
        break;
    case 2: // with Cell
        stepWriter(s);
        file << '\n';
        for(const auto& vec: s.getCellVec()*s.getCellDim(AtomFmt::Angstrom)){
            file << vec[0] << ' ' << vec[1] << ' ' << vec[2] << '\n';
        }
        break;
    default:
        throw IO::Error{"XYZ: invalid filemode setting"};
    }
    return true;
}

const IO::Plugin IO::XYZ =
{
    "xyz",
    "xyz",
    "xyz",
    &XYZParser,
    &XYZWriter,
    nullptr,
    &makePreset
};
