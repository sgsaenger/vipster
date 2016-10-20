#include "xyz.h"
#include <sstream>
#include <iostream>
#include <cstring>

using namespace Vipster;

std::tuple<Molecule,optional<Param>> Vipster::IO::XYZ_parser(std::string fn, std::ifstream &file)
{
    Molecule m{fn,1};
    char linebuf[BUFFLEN];
    char *token;
    std::vector<std::string> tokens;
    int nat;
    std::string name;
    if(!file.getline(linebuf,BUFFLEN)||!std::sscanf(linebuf,"%d",&nat)){
        throw IOError("XYZ: Could not read number of atoms.");
    }
    m.curStep().newAtoms(nat);
    if(!file.getline(linebuf,BUFFLEN)){
        throw IOError("XYZ: Could not read comment.");
    }
    m.curStep().comment = std::string(linebuf);
    for(int i=0;i!=nat;++i){
        if(!file.getline(linebuf,BUFFLEN)){
            throw IOError("XYZ: Should be "+std::to_string(nat)+" atoms, but found only:"+std::to_string(i));
        }
        token = strtok(linebuf," ");
        while(token){
            tokens.emplace_back(token);
            token = strtok(NULL," ");
        }
        m.curStep().setAtom(i,tokens[0],{std::stof(tokens[1]),std::stof(tokens[2]),std::stof(tokens[3])});
        tokens.clear();
    }

    return std::tuple<Molecule,optional<Param>>(m,optional<Param>());
}

void    Vipster::IO::XYZ_writer(const Molecule &m, std::ofstream &file,Param &p)
{

}
