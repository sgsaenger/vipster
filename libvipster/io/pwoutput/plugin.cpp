#include "plugin.h"

#include <sstream>

using namespace Vipster;

IO::Data PWOutParser(const std::string& name, std::ifstream &file)
{
    IO::Data d{};
    Molecule& m = d.mol;
    m.setName(name);
    Step *s = &m.newStep();

    std::string line, dummy_s;
    size_t nat{}, ntype{};
    float celldim{};
    Mat cellvec;
    bool gamma{false}, readInitial{false};
    CdmFmt cdmfmt{CdmFmt::Bohr};
    while (std::getline(file, line)) {
        if (line.find("number of atoms/cell") != std::string::npos) {
            std::stringstream{line.substr(33)} >> nat;
            s->newAtoms(nat);
            std::getline(file, line);
            std::stringstream{line.substr(33)} >> ntype;
        } else if (line.find("gamma-point") != std::string::npos) {
            gamma = true;
        } else if (line.find("number of k points=") != std::string::npos) {
            if (gamma) {
                continue;
            }
            KPoints kpts{};
            size_t nk = static_cast<size_t>(std::stoi(line.substr(line.find('=')+1)));
            // skip header
            std::getline(file, line);
            if (line.find("cart. coord.") == std::string::npos){
                continue;
            }
            kpts.discrete.kpoints.resize(nk);
            for (auto& k: kpts.discrete.kpoints) {
                std::getline(file, line);
                std::stringstream ss{line};
                std::getline(ss, dummy_s, '(');
                std::getline(ss, dummy_s, '(');
                ss >> k.pos[0] >> k.pos[1] >> k.pos[2];
                std::getline(ss, dummy_s, '=');
                ss >> k.weight;
            }
        } else if (line.find("celldm(1)=") != std::string::npos) {
            /*
             * parse celldm(1)
             * always given with discrete vectors
             * ibrav and rest of celldm not needed
             */
            std::stringstream{line} >> dummy_s >> celldim;
            s->setCellDim(celldim, cdmfmt);
            // skip to cell-vectors
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            for(size_t i=0; i<3; ++i){
                std::getline(file, line);
                std::stringstream{line} >> dummy_s >> dummy_s >> dummy_s
                        >> cellvec[i][0] >> cellvec[i][1] >> cellvec[i][2];
            }
            s->setCellVec(cellvec);
        } else if (!readInitial && (line.find("site n.") != std::string::npos)) {
            // parse initial coordinates
            // always given as ALAT (or aditionally as CRYSTAL with high verbosity)
            s->setFmt(AtomFmt::Alat);
            for(auto& at: *s){
                std::getline(file, line);
                std::stringstream ss{line};
                ss >> dummy_s >> at.name >> dummy_s;
                std::getline(ss, dummy_s, '(');
                ss >> at.coord[0] >> at.coord[1] >> at.coord[2];
                if (ss.fail()) {
                    throw IO::Error{"Failed to parse atom"};
                }
            }
            readInitial = true;
        } else if ((line.find("CELL_PARAMETERS") != std::string::npos) &&
                   (line.find("DEPRECATED") == std::string::npos)) {
            if (line.find("(bohr)") != std::string::npos) {
                cdmfmt = CdmFmt::Bohr;
                celldim = 1;
            }else if (line.find("angstrom") != std::string::npos) {
                cdmfmt = CdmFmt::Angstrom;
                celldim = 1;
            }else{
                cdmfmt = CdmFmt::Bohr;
                celldim = std::stof(line.substr(line.find('=')+1));
            }
            // parse vectors
            for(Vec& v: cellvec){
                std::getline(file, line);
                std::stringstream{line} >> v[0] >> v[1] >> v[2];
            }
        } else if (line.find("ATOMIC_POSITIONS") != std::string::npos) {
            // formatted positions
            // creating new step here
            s = &m.newStep();
            s->setCellDim(celldim, cdmfmt);
            s->setCellVec(cellvec);
            s->newAtoms(nat);
            if(line.find("angstrom") != std::string::npos) {
                s->setFmt(AtomFmt::Angstrom);
            }else if (line.find("bohr") != std::string::npos) {
                s->setFmt(AtomFmt::Bohr);
            }else if (line.find("crystal") != std::string::npos) {
                s->setFmt(AtomFmt::Crystal);
            }else{
                s->setFmt(AtomFmt::Alat);
            }
            for(auto& at: *s){
                std::getline(file, line);
                std::stringstream{line} >> at.name
                        >> at.coord[0] >> at.coord[1] >> at.coord[2];
            }
        }
        else if (line.find("Begin final coordinates") != std::string::npos) {
            break;
        }
    }

    return d;
}

const IO::Plugin IO::PWOutput =
{
    "PWScf Output File",
    "pwo",
    "pwo",
    IO::Plugin::None,
    &PWOutParser
};
