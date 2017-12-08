#include "ioplugins/lmpinput.h"

//TODO: REWORK, WRITER
#include <cstring>

using namespace Vipster;

std::shared_ptr<IO::BaseData> LmpInpParser(std::string name, std::ifstream &file)
{
    enum class ParseMode{Header,Atoms,Types};

    auto data = std::make_shared<IO::BaseData>();
    Molecule& m = data->mol;
    m.setName(name);
    Step& s = m.newStep();

    char line[IO::linelen];
    char* tok;
    size_t nat, ntyp, count;
    Mat cell;
    std::vector<int> typemap;
    std::vector<std::string> types;
    float t1,t2;
    auto mode = ParseMode::Header;
    while(file.getline(line, IO::linelen)){
        if(mode == ParseMode::Header){
            if(strstr(line, "atoms")){
                sscanf(line, "%ld", &nat);
                s.newAtoms(nat);
                typemap.resize(nat);
            }else if(strstr(line, "atom types")){
                sscanf(line, "%ld", &ntyp);
                types.reserve(ntyp);
            }else if(strstr(line, "xlo xhi")){
                sscanf(line, "%g %g", &t1, &t2);
                cell[0][0] = t2 - t1;
            }else if(strstr(line, "ylo yhi")){
                sscanf(line, "%g %g", &t1, &t2);
                cell[1][1] = t2 - t1;
            }else if(strstr(line, "zlo zhi")){
                sscanf(line, "%g %g", &t1, &t2);
                cell[2][2] = t2 - t1;
            }else if(strstr(line, "xy xz yz")){
                sscanf(line, "%g %g %g", &cell[0][1], &cell[0][2], &cell[1][2]);
            }else if(strstr(line, "Atoms")){
                mode = ParseMode::Atoms;
                count = -1;
            }else if(strstr(line, "Masses")){
                mode = ParseMode::Types;
                count = -1;
            }
        }else if(mode == ParseMode::Atoms){
            if(!count) file.getline(line, IO::linelen);
            Atom& at = s.getAtomMod(count-1);
            const char* fullfmt = "%*d %*d %d %f %f %f %f\n";
            sscanf(line, fullfmt, &typemap[count-1], &at.charge, &at.coord[0], &at.coord[1], &at.coord[2]);
            count++;
            if(count==nat){
                mode = ParseMode::Header;
            }
        }else if(mode == ParseMode::Types){
            if(!count){
                types.push_back("X");
                file.getline(line, IO::linelen);
            }
            tok = strstr(line, "#");
            if(tok){
                tok++;
                while(isblank(tok[0])){
                    tok++;
                }
                types.push_back(std::string(tok));
            }else{
                throw IOError("Lammps data files must specify atom types in comment to masses :(");
            }
            count++;
            if(count==ntyp){
                mode = ParseMode::Header;
            }
        }
    }
    for(size_t i = 0; i < s.getNat(); ++i){
        s.getAtomMod(i).name = types[typemap[i]];
    }
    s.setFmt(AtomFmt::Angstrom);
    s.setCellDim(1);
    s.setCellVec(cell);
    return data;
}

bool LmpInpWriter(const Molecule& m, std::ofstream &file, const IO::BaseParam* p)
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
