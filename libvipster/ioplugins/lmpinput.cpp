#include "ioplugins/lmpinput.h"

//TODO: REWORK, WRITER
#include <sstream>

using namespace Vipster;

std::shared_ptr<IO::BaseData> LmpInpParser(std::string name, std::ifstream &file)
{
    enum class ParseMode{Header,Atoms,Types};

    auto data = std::make_shared<IO::BaseData>();
    Molecule& m = data->mol;
    m.setName(name);
    StepProper& s = m.newStep();
    s.setFmt(AtomFmt::Angstrom);

    std::string line;
    size_t nat, ntype;
    float t1, t2;
    Mat cell;
    std::map<std::string, std::string> types;
    bool hasNames{false};
    while (std::getline(file, line)) {
        if (line.find("atoms") != line.npos) {
            std::stringstream ss{line};
            ss >> nat;
            if (ss.fail()) {
                throw IOError("Lammps Input: failed to"
                              "parse number of atoms");
            } else {
                s.newAtoms(nat);
            }
        } else if (line.find("atom types") != line.npos) {
            std::stringstream ss{line};
            ss >> ntype;
            if (ss.fail()) {
                throw IOError("Lammps Input: failed to"
                              "parse number of types");
            }
        } else if (line.find("xlo xhi") != line.npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IOError("Lammps Input: failed to"
                              "parse cell X dimension");
            } else {
                cell[0][0] = t2 - t1;
            }
        } else if (line.find("ylo yhi") != line.npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IOError("Lammps Input: failed to"
                              "parse cell Y dimension");
            } else {
                cell[1][1] = t2 - t1;
            }
        } else if (line.find("zlo zhi") != line.npos) {
            std::stringstream ss{line};
            ss >> t1 >> t2;
            if (ss.fail()) {
                throw IOError("Lammps Input: failed to"
                              "parse cell Z dimension");
            } else {
                cell[2][2] = t2 - t1;
            }
        } else if (line.find("xy xz yz") != line.npos) {
            std::stringstream ss{line};
            ss >> cell[1][0] >> cell[2][0] >> cell[2][1];
            if (ss.fail()) {
                throw IOError("Lammps Input: failed to"
                              "parse cell tilt factors");
            }
        } else if (line.find("Masses") != line.npos) {
            std::getline(file, line);
            std::string id, name;
            for (size_t i=0; i<ntype; ++i) {
                std::getline(file, line);
                std::stringstream ss{line};
                if (line.find("#") != line.npos) {
                    // if there's a comment, extract the typename from it
                    hasNames = true;
                    ss >> id >> t1 >> name;
                    name.erase(0, 1);
                    types[id] = name;
                } else {
                    // else just number the types accordingly
                    ss >> name >> t1;
                }
                if (ss.fail())
                    throw IOError("Lammps Input: failed to parse atom type");
                (*s.pse)[name].m = t1;
            }
        } else if (line.find("Atoms") != line.npos) {
            //TODO: format-support
            std::getline(file, line);
            // full-fmt parser only atm
            auto parse = [](std::string& s, AtomRef at){
                int d1, d2;
                std::stringstream ss{s};
                ss >> d1 >> d2 >> at.name >> at.charge
                   >> at.coord[0] >> at.coord[1] >> at.coord[2];
                if (ss.fail())
                    throw IOError("Lammps Input: failed to parse atom");
            };
            for (size_t i=0; i<nat; ++i) {
                std::getline(file, line);
                parse(line, s[i]);
            }
        }
    }
    if (hasNames) {
        for (size_t i=0; i<nat; ++i) {
            auto at = s[i];
            at.name = types[at.name];
        }
    }
    return data;
}

bool LmpInpWriter(const Molecule&, std::ofstream &, const IO::BaseParam* p)
{
    auto *lp = dynamic_cast<const IO::LmpParam*>(p);
    if(!lp) throw IOError("Lammps-Writer needs parameter set");
    return false;
}

const IOPlugin IO::LmpInput =
{
    "Lammps Data File",
    "lmp",
    "lmp",
    &LmpInpParser,
    nullptr
};
