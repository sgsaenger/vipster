#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>
#include "molecule.h"
#include "ioplugin.h"
#include "ioplugins/xyz.h"
#include "ioplugins/pwinput.h"
#include "ioplugins/pwoutput.h"
#include "ioplugins/lmpinput.h"
#include "ioplugins/lmptrajec.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    const std::map<IOFmt, IO::Plugin const *const> IOPlugins{
            {IOFmt::XYZ, &IO::XYZ},
            {IOFmt::PWI, &IO::PWInput},
            {IOFmt::PWO, &IO::PWOutput},
            {IOFmt::LMP, &IO::LmpInput},
            {IOFmt::DMP, &IO::LmpTrajec}
    };
    IO::Data readFile(const std::string &fn, IOFmt fmt);
    IO::Data readFile(const std::string &fn, IOFmt fmt, std::string name);
    bool     writeFile(const std::string &fn, IOFmt fmt, const Molecule &m,
                       const BaseParam *p=nullptr,
                       const BaseConfig *c=nullptr);
    /*
     * Convenience maps
     */
    const auto IOCmdIn = [](){
        std::map<std::string, IOFmt> fmts_in;
        for(const auto& pair: IOPlugins){
            fmts_in[pair.second->command] = pair.first;
        }
        return fmts_in;
    }();
    const auto IOCmdOut = [](){
        std::map<std::string, IOFmt> fmts_out;
        for(const auto& pair: IOPlugins){
            if(pair.second->writer!=nullptr){
                fmts_out[pair.second->command] = pair.first;
            }
        }
        return fmts_out;
    }();
    const auto IOExt = [](){
        std::map<std::string, IOFmt> ext;
        for(const auto& pair: IOPlugins){
            ext[pair.second->extension] = pair.first;
        }
        return ext;
    }();
}

#endif // IOWRAPPER

