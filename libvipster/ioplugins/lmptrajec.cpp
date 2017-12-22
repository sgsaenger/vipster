#include "ioplugins/lmptrajec.h"

//TODO: REWORK
#include <cstring>

using namespace Vipster;
using parseFunc = void (*)(const char*, Atom& at);

parseFunc IdentifyColumn(const char *const tok){
    constexpr char const* names[] = {
        "element",
        "x", "y", "z",
//        "xs", "ys", "zs",
        "q"
    };
    const parseFunc funcMap[] = {
        [](const char* c, Atom& at){at.name = c;},
        [](const char* c, Atom& at){at.coord[0] = atof(c);},
        [](const char* c, Atom& at){at.coord[1] = atof(c);},
        [](const char* c, Atom& at){at.coord[2] = atof(c);},
        [](const char* c, Atom& at){at.charge = atof(c);},
    };
    for(size_t i=0; i<5; ++i){
        if(tok[0] == names[i][0]) return funcMap[i];
    }
    return [](const char*, Atom&){};
}

std::shared_ptr<IO::BaseData> LmpTrajecParser(std::string name, std::ifstream &file)
{
    enum class ParseMode{Header, Cell, Atoms};
    constexpr int xl[]={0,0,1}, yl[]={1,2,2};

    auto data = std::make_shared<IO::BaseData>();
    Molecule& m = data->mol;
    m.setName(name);
    StepProper* s = nullptr;

    char line[IO::linelen];
    Mat cell;
    float t1, t2, t3;
    int nat, count;
    bool firstpass = true;
    ParseMode mode = ParseMode::Header;
    auto atomfmt = AtomFmt::Angstrom;
    std::vector<parseFunc> parseTable;
    while(file.getline(line, IO::linelen)){
        if(mode == ParseMode::Header){
            if(strstr(line, "TIMESTEP")){
                file.getline(line, IO::linelen);
                s = &m.newStep();
                s->setFmt(atomfmt);
            }else if(strstr(line, "NUM")){
                file.getline(line, IO::linelen);
                sscanf(line, "%d", &nat);
                s->newAtoms(nat);
            }else if(strstr(line, "BOX")){
                mode = ParseMode::Cell;
                count = 0;
                cell = Mat{};
            }else if(strstr(line, "ATOMS")){
                mode = ParseMode::Atoms;
                count = 0;
                if (firstpass) {
                    if(strstr(line, "xs")){
                        atomfmt = AtomFmt::Crystal;
                        s->setFmt(atomfmt);
                    }
                    strtok(line, " ");
                    strtok(nullptr, " ");
                    char* tok;
                    while((tok = strtok(nullptr, " "))){
                        parseTable.push_back(IdentifyColumn(tok));
                    }
                    firstpass = false;
                }
            }
        }else if(mode == ParseMode::Cell){
            int test = sscanf(line, "%f %f %f", &t1, &t2, &t3);
            cell[count][count] = t2 - t1;
            if (test == 3) cell[xl[count]][yl[count]] = t3;
            count++;
            if (count == 3){
                s->setCellVec(cell);
                mode = ParseMode::Header;
            }
        }else if(mode == ParseMode::Atoms){
            auto at = (*s)[count];
            parseTable[0](strtok(line, " "), at);
            int i = 0;
            char *tok;
            while((tok = strtok(nullptr," "))){
                parseTable[++i](tok, at);
            }
            count++;
            if (count == nat){
                mode = ParseMode::Header;
            }
        }
    }
    return data;
}

const IOPlugin IO::LmpTrajec =
{
    "Lammps Dump",
    "lammpstrj",
    "dmp",
    &LmpTrajecParser,
    nullptr
};
