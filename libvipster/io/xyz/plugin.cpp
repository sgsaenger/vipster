#include "plugin.h"

#include <sstream>
#include <iomanip>

using namespace Vipster;

IO::Data XYZParser(const std::string& name, std::istream &file)
{
    IO::Data data{};
    Molecule &m = data.mol;
    m.setName(name);

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
                        m.getStep(0).setCellDim(1, CdmFmt::Angstrom);
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
               const IO::BaseParam*const, const IO::BasePreset*const c,
               size_t index)
{
    const auto *cc = dynamic_cast<const IO::XYZPreset*>(c);
    if(!cc) throw IO::Error("XYZ-Writer needs IO preset");
    const Step& s = m.getStep(index).asFmt(AtomFmt::Angstrom);
    auto stepWriter = [&file, cc](const Step& s){
        file << s.getNat() << '\n';
        file << s.getComment() << '\n';
        file << std::fixed << std::setprecision(5);
        switch(cc->atomdata){
        case IO::XYZPreset::Data::None:
            for(const auto& at: s){
                file << std::left << std::setw(3) << at.name << " "
                     << std::right << std::setw(10) << at.coord[0] << " "
                     << std::right << std::setw(10) << at.coord[1] << " "
                     << std::right << std::setw(10) << at.coord[2] << '\n';
            }
            break;
        case IO::XYZPreset::Data::Charge:
            for(const auto& at: s){
                file << std::left << std::setw(3) << at.name << " "
                     << std::right << std::setw(10) << at.coord[0] << " "
                     << std::right << std::setw(10) << at.coord[1] << " "
                     << std::right << std::setw(10) << at.coord[2] << " "
                     << std::right << std::setw(10) << at.properties->charge << '\n';
            }
            break;
        case IO::XYZPreset::Data::Forces:
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
        }
    };
    switch(cc->filemode){
    case IO::XYZPreset::Mode::Step:
        stepWriter(s);
        break;
    case IO::XYZPreset::Mode::Cell:
        stepWriter(s);
        file << '\n';
        for(const auto& vec: s.getCellVec()*s.getCellDim(CdmFmt::Angstrom)){
            file << vec[0] << ' ' << vec[1] << ' ' << vec[2] << '\n';
        }
        break;
    case IO::XYZPreset::Mode::Trajec:
        for(const auto& step: m.getSteps()){
            stepWriter(step.asFmt(AtomFmt::Angstrom));
        }
        break;
    }
    return true;
}

static std::unique_ptr<IO::BasePreset> makePreset()
{
    return std::make_unique<IO::XYZPreset>();
}

const IO::Plugin IO::XYZ =
{
    "xyz",
    "xyz",
    "xyz",
    IO::Plugin::Read | IO::Plugin::Write | IO::Plugin::Preset,
    &XYZParser,
    &XYZWriter,
    nullptr,
    &makePreset
};
