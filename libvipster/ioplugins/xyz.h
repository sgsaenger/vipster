#ifndef XYZ_H
#define XYZ_H

#include "ioplugin.h"
#include <cstring>

namespace Vipster{
namespace IO{
const IOPlugin XYZ{"xyz","xyz","xyz",
    // Parser:
    [](std::string fn, std::ifstream &file) -> std::tuple<Molecule, IOType, IOBase*>{
        Molecule m{fn,0};
        char linebuf[BUFFLEN];
        char *token;
        std::vector<std::string> tokens;
        int nat;
        std::string name;
        if(!file.getline(linebuf,BUFFLEN)||!std::sscanf(linebuf,"%d",&nat)){
            throw IOError("XYZ: Could not read number of atoms.");
        }
        Step s;
        s.newAtoms(nat);
        if(!file.getline(linebuf,BUFFLEN)){
            throw IOError("XYZ: Could not read comment.");
        }
        s.comment = std::string(linebuf);
        for(int i=0;i!=nat;++i){
            if(!file.getline(linebuf,BUFFLEN)){
                throw IOError("XYZ: Should be "+std::to_string(nat)+" atoms, but found only:"+std::to_string(i));
            }
            token = strtok(linebuf," ");
            while(token){
                tokens.emplace_back(token);
                token = strtok(NULL," ");
            }
            s.setAtom(i,tokens[0],{std::stof(tokens[1]),std::stof(tokens[2]),std::stof(tokens[3])});
            tokens.clear();
        }
        m.insertStep(s);
        return std::tuple<Molecule, IOType, IOBase*>(m, IOType::None, nullptr);
    },
    // Writer:
    [](const Molecule &m, std::ofstream &file, IOBase &p) -> void{

    }
};
}
}

#endif // XYZ_H
