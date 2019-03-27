#include "plugin.h"
#include "../util.h"
#include <sstream>

using namespace Vipster;

template<typename T>
auto mkparser(){
    return [](IO::Data& d, std::string&& n, std::ifstream& f){
        constexpr size_t N = T::Dim;
        std::array<size_t, N> extent;
        f >> extent[0] >> extent[1] >> extent[2];
        T grid{extent};
        grid.name = n;
        f >> grid.origin[0] >> grid.origin[1] >> grid.origin[2];
        grid.origin *= invbohr;
        for(size_t i=0; i<N; ++i){
            f >> grid.cell[i][0] >> grid.cell[i][1] >> grid.cell[i][2];
        }
        grid.cell *= invbohr;
        grid.origin += Vec{0.5f/extent[0], 0.5f/extent[1], 0.5f/extent[2]}*grid.cell;
        for(auto& val: grid){
            f >> val;
        }
        d.data.emplace_back(std::make_unique<const T>(std::move(grid)));
    };
};

IO::Data XSFParser(const std::string& name, std::ifstream &file)
{
    IO::Data data{};
    auto& m = data.mol;
    m.setName(name);

    std::string line, buf;
    bool animated{false};
    bool hasCell{false};
    std::ifstream::pos_type fpos;
    size_t dim{0}, grp{1};
    auto tokenizeLine = [&](){
        std::getline(file, line);
        std::stringstream ls{line};
        std::vector<std::string> tokens{std::istream_iterator<std::string>{ls}, {}};
        return tokens;
    };
    auto makeAtom = [](Step& s, const std::vector<std::string>& tokens){
        AtomProperties prop{};
        if(tokens.size() == 7){
            prop.forces = Vec{stof(tokens[4]),
                              stof(tokens[5]),
                              stof(tokens[6])};
        }
        s.newAtom(tokens[0],
                  Vec{stof(tokens[1]),
                      stof(tokens[2]),
                      stof(tokens[3])},
                  prop);
    };

    while(std::getline(file, line)){
        if(line[0] == '#'){
            continue;
        }
        line = IO::trim(line);
        if(line.find("ANIMSTEPS") != line.npos){
            size_t nstep;
            std::stringstream ls{line};
            ls >> buf >> nstep;
            if(ls.fail()){
                throw IO::Error("XSF: Could not parse number of ANIMSTEPS");
            }
            m = Molecule{name, nstep};
            animated = true;
        }else if(line.find("CRYSTAL") != line.npos){
            dim = 3;
            hasCell = true;
        }else if(line.find("SLAB") != line.npos){
            dim = 2;
            hasCell = true;
        }else if(line.find("POLYMER") != line.npos){
            dim = 1;
            hasCell = true;
        }else if(line.find("MOLECULE") != line.npos){
            dim = 0;
            hasCell = true;
        }else if(line.find("DIM-GROUP") != line.npos){
            hasCell = true;
            std::getline(file, buf);
            std::stringstream ls{buf};
            ls >> dim >> grp;
            if(ls.fail()){
                throw IO::Error("XSF: Could not parse DIM-GROUP");
            }
        }else if(!hasCell && (line.find("ATOMS") != line.npos)){
            size_t idx{0};
            if(animated){
                size_t tmp;
                std::stringstream ls{line};
                ls >> buf >> tmp;
                if(ls.fail()){
                    throw IO::Error("XSF: Could not parse current index");
                }
                idx = tmp-1;
            }else{
                m.newStep();
            }
            auto& s = m.getStep(idx);
            s.setFmt(AtomFmt::Angstrom);
            std::vector<std::string> tokens;
            // check next line for atom, save last pos to rewind when going too far
            while((void)(fpos = file.tellg()),
                  (void)(tokens = tokenizeLine()),
                  (tokens.size() == 4)  || (tokens.size() == 7)){
                makeAtom(s, tokens);
            }
            // rewind last line
            file.seekg(fpos);
        }else if(hasCell && (line.find("PRIMCOORD") != line.npos)){
            size_t idx{0};
            if(animated){
                size_t tmp;
                std::stringstream ls{line};
                ls >> buf >> tmp;
                if(ls.fail()){
                    throw IO::Error("XSF: Could not parse current index");
                }
                idx = tmp-1;
            }
            auto& s = m.getStep(idx);
            s.setFmt(AtomFmt::Angstrom);
            std::getline(file, line);
            size_t nat{std::stoul(line)};
            for(size_t i=0; i<nat; ++i){
                makeAtom(s, tokenizeLine());
            }
        }else if(hasCell && (line.find("VEC") != line.npos)){
            Mat tmp;
            for(auto& v: tmp){
                file >> v[0] >> v[1] >> v[2];
            }
            if(hasCell && (line.find("CONVVEC") != line.npos)){
                // ignore conventional cell for now
                // would be useful, but is optional
                continue;
            }
            if(animated){
                size_t idx;
                std::stringstream ls{line};
                ls >> buf >> idx;
                if(ls.fail()){
                    // apply to whole trajectory
                    for(auto& s: m.getSteps()){
                        s.setCellVec(tmp);
                        s.setCellDim(1, CdmFmt::Angstrom);
                    }
                }else{
                    auto& s = m.getStep(idx-1);
                    s.setCellVec(tmp);
                    s.setCellDim(1, CdmFmt::Angstrom);
                }
            }else{
                auto& s = m.newStep();
                s.setCellVec(tmp);
                s.setCellDim(1, CdmFmt::Angstrom);
            }
        }else if(hasCell && (line.find("CONVCOORD") != line.npos)){
         // ignored
            std::getline(file, line);
            size_t nat{std::stoul(line)};
            for(size_t i=0; i<nat; ++i){
                std::getline(file, line);
            }
        }else if(line.find("BEGIN_") != line.npos){
            auto type = line.substr(6);
            std::string block_end, data_begin, data_end;
            std::function<void(IO::Data&, std::string&&, std::ifstream&)> parser;
            if(type == "INFO"){
                // ignore
                block_end = "END_INFO";
                data_begin.clear();
                data_end.clear();
                parser = nullptr;
            }else if(type == "BLOCK_BANDGRID_3D"){
                // ignore for now
                block_end = "END_BLOCK_BANDGRID_3D";
                data_begin.clear();
                data_end.clear();
                parser = nullptr;
//                data_begin = "BANDGRID_3D";
//                data_end = "END_BANDGRID_3D";
            }else if(type == "BLOCK_DATAGRID_2D"){
                // make DataGrid2_f
                block_end = "END_BLOCK_DATAGRID_2D";
                data_begin = "DATAGRID_2D";
                data_end = "END_DATAGRID_2D";
                parser = mkparser<DataGrid2D_f>();
            }else if(type == "BLOCK_DATAGRID_3D"){
                // make DataGrid3_f
                block_end = "END_BLOCK_DATAGRID_3D";
                data_begin = "DATAGRID_3D";
                data_end = "END_DATAGRID_3D";
                parser = mkparser<DataGrid3D_f>();
            }else{
                throw IO::Error("XSF: unknown block: " + line);
            }
            if(!data_begin.empty()){
                while((void)(file >> buf), buf!=block_end){
                    if(buf.find(data_begin) != buf.npos){
                        parser(data, name + ' ' + buf.substr(data_begin.size()+1), file);
                        file >> buf;
                        if(buf.find(data_end) == buf.npos){
                            throw IO::Error("XSF: Invalid datagrid" + buf);
                        }
                    }
                }
            }
        }
    }

    return data;
}

const IO::Plugin IO::XSF =
{
    "XCrysDen Structure File",
    "xsf",
    "xsf",
    IO::Plugin::None,
    &XSFParser
};
