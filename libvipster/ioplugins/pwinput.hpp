#ifndef PWINPUT_HPP
#define PWINPUT_HPP

#include <ioplugin.h>
#include <sstream>

namespace Vipster {
namespace IO{

typedef std::array<std::vector<std::pair<std::string, std::string>, 5> PWParam;

const IOPlugin PWInput{"PWScf Input","pwi","pwi",
    //Parser:
    [](std::string fn, std::ifstream &file) -> std::tuple<Molecule, IOType, IOBase*>{
        Molecule m{fn};
        Step& s = m.steps[0];
        enum class ParseStage{none,nl,species,positions,
                              kpoints,cell,constraints,occupations,forces};
        ParseStage stage{ParseStage::none};
        int i=0,j=0;
        int nl = 0;
        int nat = 0;
        int ntype = 0;
        std::array<Vec, 3> tvec;
        PWParam param;
        std::string line;
        while(std::getline(file, line) &&
              (line.find_first_not_of(" \n\v\f\r\t")!=std::string::npos)){
            switch(stage){
            case ParseStage::nl:
                break;
            case ParseStage::species:
                break;
            case ParseStage::positions:
                break;
            case ParseStage::kpoints:
                break;
            case ParseStage::cell:
                break;
            case ParseStage::none:
                if(line == "&CONTROL"){ stage = ParseStage::nl; nl=0;
                }else if(line == "&SYSTEM"){ stage = ParseStage::nl; nl=1;
                }else if(line == "&ELECTRONS"){ stage = ParseStage::nl; nl=2;
                }else if(line == "&IONS"){ stage = ParseStage::nl; nl=3;
                }else if(line == "&CELL"){ stage = ParseStage::nl; nl=4;
                }else if(line == "ATOMIC_SPECIES"){ stage = ParseStage::species;
                }else if(line.compare(0,16,"ATOMIC_POSITIONS")){ stage = ParseStage::positions;
                }else if(line.compare(0,8,"K_POINTS")){ stage = ParseStage::kpoints;
                }else if(line.compare(0,15,"CELL_PARAMETERS")){ stage = ParseStage::cell;
                }else if(line == "OCCUPATIONS"){ stage = ParseStage::occupations;
                }else if(line == "CONSTRAINTS"){ stage = ParseStage::constraints;
                }else if(line == "ATOMIC_FORCES"){ stage = ParseStage::forces;}
                break;
            }
        }
    },
    //Writer:
    [](const Molecule &m, std::ofstream &file, IOBase &p) -> void{

    }
};
}
}

#endif // PWINPUT_HPP
