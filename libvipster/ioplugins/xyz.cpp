#include "ioplugins/xyz.h"

#include <iomanip>

using namespace Vipster;

std::shared_ptr<IO::BaseData> XYZParser(std::string name, std::ifstream &file)
{
    auto data = std::make_shared<IO::BaseData>();
    Molecule &m = data->mol;
    m.setName(name);

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream natline{line};
        int nat;
        natline >> nat;
        if (natline.fail()) {
            if (!m.getNstep()) throw IOError("XYZ: Failed to parse nat");
            else {
                while ((natline>>nat).fail()) {
                    if (file.eof())
                        throw IOError("XYZ: Non-standard data after XYZ-file");
                    natline = std::stringstream{line};
                }
            }
        }
        StepProper &sp = m.newStep();
        sp.setFmt(AtomFmt::Angstrom);
        sp.newAtoms(nat);
        std::getline(file, line);
        sp.setComment(line);
        for (auto& at: sp) {
            std::getline(file, line);
            std::stringstream atline{line};
            atline >> at.name >> at.coord[0] >> at.coord[1] >> at.coord[2];
            if(atline.fail()) throw IOError("XYZ: failed to parse atom");
        }
    }
    return data;
}

bool XYZWriter(const Molecule& m, std::ofstream &file, const IO::BaseParam*)
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

const IOPlugin IO::XYZ =
{
    "xyz",
    "xyz",
    "xyz",
    &XYZParser,
    &XYZWriter
};
