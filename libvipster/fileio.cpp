#include <iostream>
#include "fileio.h"
#include "io/xyz/plugin.h"
#include "io/pwinput/plugin.h"
#include "io/pwoutput/plugin.h"
#include "io/lmpinput/plugin.h"
#include "io/lmptrajec/plugin.h"
#include "io/cpmdinput/plugin.h"
#include "io/cube/plugin.h"
#include "io/xsf/plugin.h"
#include "io/orca/plugin.h"
#include "io/poscar/plugin.h"

using namespace Vipster;

IO::Plugins IO::defaultPlugins()
{
    return {
            &IO::XYZ,
            &IO::PWInput,
            &IO::PWOutput,
            &IO::LmpInput,
            &IO::LmpTrajec,
            &IO::CPInput,
            &IO::Cube,
            &IO::XSF,
            &IO::OrcaInput,
            &IO::Poscar
    };
}

IO::Parameters Vipster::IO::defaultParams(const Plugins &p)
{
    IO::Parameters tmp;
    for(const auto& plug: p){
        if(plug->makeParam){
            tmp[plug]["default"] = plug->makeParam();
        }
    }
    return tmp;
}

IO::Presets Vipster::IO::defaultPresets(const Plugins &p)
{
    IO::Presets tmp;
    for(const auto& plug: p){
        if(plug->makePreset){
            tmp[plug]["default"] = plug->makePreset();
        }
    }
    return tmp;
}

std::optional<const IO::Plugin*> Vipster::guessFmt(std::string fn,
                                                   const IO::Plugins &p)
{
    auto pos = fn.find_last_of('.');
    if(pos != fn.npos){
        fn = fn.substr(pos+1);
    }else{
        pos = fn.find_last_of("/\\");
        if(pos != fn.npos){
            fn = fn.substr(pos+1);
        }
    }
    auto plug = std::find_if(p.begin(), p.end(), [&](const IO::Plugin* p){
        return p->extension == fn;
    });
    if(plug != p.end()){
        return *plug;
    }else{
        return {};
    }
}

// read with format guess
IO::Data Vipster::readFile(const std::string &fn, const IO::Plugins &p)
{
    // get format
    auto plugin = guessFmt(fn, p);
    if(!plugin){
        throw IO::Error{"Could not deduce format of file "+fn+
                        "\nPlease specify format explicitely", false};
    }
    // read file
    return readFile(fn, plugin.value());
}

// read with explicit format
IO::Data Vipster::readFile(const std::string &fn, const IO::Plugin *plug)
{
    // check if file can be read
    std::ifstream file{fn};
    if(!file){
        throw IO::Error("Could not open "+fn);
    }
    // try to parse
    if(!plug->parser){
        throw IO::Error("Format is not readable");
    }
    auto tmp = plug->parser(fn, file);
    if(!tmp.mol.getNstep()){
        throw IO::Error("No Molecule could be parsed");
    }
    // return if successful
    return tmp;
}

bool  Vipster::writeFile(const std::string &fn,
                         const IO::Plugin *plug,
                         const Molecule &m,
                         std::optional<size_t> idx,
                         const IO::BaseParam *const p,
                         const std::optional<IO::BasePreset>& c)
{
    std::ofstream file{fn};
    if(!idx){
        idx = m.getNstep()-1;
    }
    try{
        if(!file){
            throw IO::Error{"Could not open "+fn};
        }
        if(!plug->writer){
            throw IO::Error{"Read-only format"};
        }
        return plug->writer(m, file, p, c, *idx);
    }
    catch(IO::Error &e){
        std::cout << e.what() << std::endl;
        throw;
    }
}
