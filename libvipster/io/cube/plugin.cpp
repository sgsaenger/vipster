#include "plugin.h"

using namespace Vipster;

IO::Data CubeParser(const std::string& name, std::ifstream &file){
    IO::Data d{};
    d.fmt = IOFmt::CUBE;
    d.mol.setName(name);
    StepProper &s = d.mol.newStep();

    std::string line, buf;
    std::stringstream lstream;

    // two comment lines
    std::getline(file, buf);
    std::getline(file, line);
    if(!buf.empty() && !line.empty()){
        s.setComment(buf + " | " + line);
    }else if(!buf.empty()){
        s.setComment(buf);
    }else if(!line.empty()){
        s.setComment(line);
    }

    std::getline(file, line);
    lstream.str(line);
    long tmp;
    bool mo_mode{false};
    size_t nat, nval;
    Vec origin;
    Mat cell;
    std::array<size_t, 3> extent;
    lstream >> tmp >> origin[0] >> origin[1] >> origin[2];
    if(lstream.fail()){
        throw IO::Error("Cube: failed to parse nat/origin");
    }
    nat = static_cast<size_t>(std::abs(tmp));
    if(tmp<0){
        mo_mode = true;
    }
    s.newAtoms(static_cast<size_t>(nat));
    lstream >> nval;
    if(lstream.fail()){
        nval = 1;
    }
    for(size_t i=0; i<3; ++i){
        file >> tmp >> cell[i][0] >> cell[i][1] >> cell[i][2];
        extent[i] = static_cast<size_t>(std::abs(tmp));
        if(tmp<0){
            cell[i] *= extent[i]*bohrrad;
        }else{
            cell[i] *= extent[i];
        }
    }
    for(auto& at: s){
        file >> at.name >> at.properties->charge
             >> at.coord[0] >> at.coord[1] >> at.coord[2];
    }
    if(mo_mode){
        // twice to consume left-over '\n'
        std::getline(file, line);
        std::getline(file, line);
        lstream.str(line);
        lstream.clear();
        size_t nmo;
        lstream >> nmo;
        std::vector<size_t> mo_indices(nmo);
        std::vector<std::vector<float>> grids;
        std::vector<std::vector<float>::iterator> iters;
        for(size_t i=0; i<nmo; ++i){
            lstream >> mo_indices[i];
            grids.emplace_back(extent[0]*extent[1]*extent[2]);
            iters.emplace_back(grids.back().begin());
        }
        while(iters.back() != grids.back().end()){
            for(auto& it: iters){
                file >> *it++;
            }
        }
    }else{
        if(nval==1){
            std::vector<float> grid(extent[0]*extent[1]*extent[2]);
            for(auto& p: grid){
                file >> p;
            }
        }else if(nval==4){
            std::vector<float> s_grid(extent[0]*extent[1]*extent[2]);
            std::vector<std::array<float, 3>> v_grid(extent[0]*extent[1]*extent[2]);
            auto s_it = s_grid.begin();
            auto v_it = v_grid.begin();
            while(s_it != s_grid.end()){
                auto& v = *v_it;
                file >> *s_it >> v[0] >> v[1] >> v[2];
                s_it++;
                v_it++;
            }
        }else{
            throw IO::Error("Cube: Only scalar or scalar+gradient grids are supported");
        }
    }

    return d;
}

const IO::Plugin IO::Cube =
{
    "Gaussian Cube file",
    "cube",
    "cube",
    IO::Plugin::None,
    &CubeParser,
    nullptr,
    nullptr,
    nullptr
};
