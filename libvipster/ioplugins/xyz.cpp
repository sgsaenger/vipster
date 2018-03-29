#include "xyz.h"

#include <sstream>
#include <iomanip>

using namespace Vipster;

IO::Data XYZParser(std::string name, std::ifstream &file)
{
    IO::Data data{};
    data.fmt = IOFmt::XYZ;
    Molecule &m = data.mol;
    m.setName(name);

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream natline{line};
        size_t nat;
        natline >> nat;
        if (natline.fail()) {
            if (!m.getNstep()) throw IO::Error("XYZ: Failed to parse nat");
            else {
                while ((natline>>nat).fail()) {
                    std::getline(file, line);
                    if (file.eof())
                        throw IO::Error("XYZ: Non-standard data after XYZ-file");
                    natline = std::stringstream{line};
                }
            }
        }
        StepProper &sp = m.newStep();
        sp.setFmt(AtomFmt::Angstrom);
        sp.enableCell(false);
        sp.newAtoms(nat);
        std::getline(file, line);
        sp.setComment(line);
        for (auto& at: sp) {
            std::getline(file, line);
            std::stringstream atline{line};
            atline >> at.name >> at.coord[0] >> at.coord[1] >> at.coord[2];
            if(atline.fail()) throw IO::Error("XYZ: failed to parse atom");
        }
    }
    return data;
}

bool XYZWriter(const Molecule& m, std::ofstream &file,
               const BaseParam*const, const BaseConfig*const)
{
    const Step& s = m.getStep(0).asAngstrom;
    file << s.getNat() << '\n';
    file << s.getComment() << '\n';
    file << std::fixed << std::setprecision(5);
    for(auto& at: s){
        file << std::left << std::setw(3) << at.name << " "
             << std::right << std::setw(10) << at.coord[0] << " "
             << std::right << std::setw(10) << at.coord[1] << " "
             << std::right << std::setw(10) << at.coord[2] << '\n';
    }
    return true;
}

const IO::Plugin IO::XYZ =
{
    "xyz",
    "xyz",
    "xyz",
    IO::Plugin::None,
    &XYZParser,
    &XYZWriter
};
