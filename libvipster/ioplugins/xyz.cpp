#include "ioplugins/xyz.h"

//TODO: REWORK C-STRINGS, WRITER
#include <iostream>
#include <iomanip>

using namespace Vipster;

std::shared_ptr<IO::BaseData> XYZParser(std::string name, std::ifstream &file)
{
    enum class ParseMode{Header, Atoms};

    auto data = std::make_shared<IO::BaseData>();
    Molecule &m = data->mol;
    m.setName(name);

    auto mode = ParseMode::Header;
    int nat, count;
    char line[IO::linelen], type[10];
    StepProper *sp = nullptr;
    while(file.getline(line, IO::linelen)){
        if(mode == ParseMode::Header){
            int test = sscanf(line, "%d", &nat);
            if(test != 1) continue;
            sp = &m.newStep();
            sp->setFmt(AtomFmt::Angstrom);
            sp->newAtoms(nat);
            count = 0;
            file.getline(line, IO::linelen);
            sp->setComment(line);
            mode = ParseMode::Atoms;
        }else if(mode == ParseMode::Atoms){
            Atom at = (*sp)[count];
            sscanf(line, "%s %f %f %f", type, &at.coord[0], &at.coord[1], &at.coord[2]);
            at.name = std::string(type);
            count++;
            if(count == nat) mode = ParseMode::Header;
        }
    }
    return data;
}

bool XYZWriter(const Molecule& m, std::ofstream &file, const IO::BaseParam*)
{
    const StepProper& s = m.getStep(0);
    file << s.getNat() << '\n';
    file << s.getComment() << '\n';
    file << std::fixed << std::setprecision(5);
//    for(auto at: s.getAtoms()){
//        file << std::left << std::setw(3) << at.name << " "
//             << std::right << std::setw(10) << at.coord[0] << " "
//             << std::right << std::setw(10) << at.coord[1] << " "
//             << std::right << std::setw(10) << at.coord[2] << '\n';
//    }
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
