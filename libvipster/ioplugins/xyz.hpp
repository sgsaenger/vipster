#ifndef XYZ_H
#define XYZ_H

#include "ioplugin.h"
#include <sstream>

namespace Vipster{
namespace IO{
const IOPlugin XYZ{"xyz","xyz","xyz",
    // Parser:
    [](std::string fn, std::ifstream &file) -> std::tuple<Molecule, IOType, IOBase*>{
        Molecule m{fn,0};
        enum class ParseStage{nat, comment, atoms};
        ParseStage stage{ParseStage::nat};
        std::string line;
        int nat = 0;
        int i = 0;
        std::shared_ptr<Step> s;
        while(std::getline(file, line)){
            switch(stage){
            case ParseStage::nat:
                try{
                    nat = std::stoi(line);
                }catch(...){
                    if(m.steps.size()){
                        break;
                    }else{
                        throw IOError("XYZ: Couldn't parse the number of atoms");
                    }
                }
                s = std::make_shared<Step>();
                stage = ParseStage::comment;
                break;
            case ParseStage::comment:
                s->comment = line;
                stage = ParseStage::atoms;
                i = 0;
                break;
            case ParseStage::atoms:
                Atom a{"C", {0.,0.,0.}, 0., {false, false, false}, false};
                std::istringstream linestream{std::move(line)};
                linestream >> a.name >> a.coord[0] >> a.coord[1] >> a.coord [2];
                s->newAtom(std::move(a));
                i++;
                if(i == nat){
                    stage = ParseStage::nat;
                    s->setCellDim(1,true,AtomFmt::Angstrom);
                    m.insertStep(std::move(*s));
                }
                break;
            }
        }
        return std::tuple<Molecule, IOType, IOBase*>(m, IOType::None, nullptr);
    },
    // Writer:
    [](const Molecule &m, std::ofstream &file, IOBase &p) -> void{

    }
};
}
}

#endif // XYZ_H
