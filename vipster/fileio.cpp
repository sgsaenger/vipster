#include <iostream>
#include <locale>
#include "fileio.h"

using namespace Vipster;
namespace fs = std::filesystem;

detail::TempWrap::TempWrap()
    : tmppath{fs::temp_directory_path()/"vipster"}
{
    if(!fs::exists(tmppath)){
        fs::create_directory(tmppath);
    }else if(!fs::is_directory(tmppath)){
        fs::remove(tmppath);
        fs::create_directory(tmppath);
    }
}

const fs::path& detail::TempWrap::getPath() const
{
    return tmppath;
}

const fs::path& Vipster::getTempPath()
{
    return detail::tempwrap.getPath();
}

const detail::TempWrap detail::tempwrap{};

const Plugin *Vipster::guessFmt(std::string fn, const PluginList &p)
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
    auto plug = std::find_if(p.begin(), p.end(), [&](const Plugin* p){
        return p->extension == fn;
    });
    if(plug != p.end()){
        return *plug;
    }else{
        return {};
    }
}

// read with explicit format
IOTuple Vipster::readFile(const std::string &fn, const Plugin *plug)
{
    if(!plug){
        throw Error{"readFile: no plugin provided"};
    }
    // set locale to C to get consistent parsing
    std::string userLocale = setlocale(0, nullptr);
    setlocale(LC_ALL, "C");
    // fail on write-only parser
    if(!plug->parser){
        setlocale(LC_ALL, userLocale.c_str());
        throw IOError("Format is not readable");
    }
    // check if file can be read
    std::ifstream file{fn};
    if(!file){
        setlocale(LC_ALL, userLocale.c_str());
        throw IOError("Could not open \""+fn+'"');
    }
    // try to parse
    auto tmp = plug->parser(fn, file);
    if(!std::get<0>(tmp).getNstep()){
        setlocale(LC_ALL, userLocale.c_str());
        throw IOError("No Molecule could be parsed");
    }
    // return if successful
    setlocale(LC_ALL, userLocale.c_str());
    return tmp;
}

IOTuple Vipster::readCin(const Plugin *plug)
{
    if(!plug){
        throw Error{"readCin: no plugin provided"};
    }
    // set locale to C to get consistent parsing
    std::string userLocale = setlocale(0, nullptr);
    setlocale(LC_ALL, "C");
    // fail on write-only parser
    if(!plug->parser){
        setlocale(LC_ALL, userLocale.c_str());
        throw IOError("Format is not readable");
    }
    // try to parse
    auto tmp = plug->parser("Molecule", std::cin);
    if(!std::get<0>(tmp).getNstep()){
        setlocale(LC_ALL, userLocale.c_str());
        throw IOError("No Molecule could be parsed");
    }
    // return if successful
    setlocale(LC_ALL, userLocale.c_str());
    return tmp;
}

bool  Vipster::writeFile(const std::string &fn,
                         const Plugin *plug,
                         const Molecule &m,
                         std::optional<size_t> idx,
                         const std::optional<Parameter>& p,
                         const std::optional<Preset>& c)
{
    if(!plug){
        throw Error{"writeFile: no plugin provided"};
    }
    if(!idx){
        idx = m.getNstep()-1;
    }
    try{
        if(!plug->writer){
            throw IOError{"Read-only format"};
        }
        bool use_temp = true;
        bool res = false;
        auto filename = getTempPath()/fs::path{fn}.filename();
        {
            // scoped writing so file is closed before copy
            std::ofstream file{filename};
            if(!file){
                use_temp = false;
                file = std::ofstream{fn};
                if(!file){
                    throw IOError{"Could not open "+fn};
                }
            }
            res = plug->writer(m, file, p, c, *idx);
        }
        if(use_temp){
            fs::copy_file(filename, fn, fs::copy_options::overwrite_existing);
        }
        return res;
    }
    catch(IOError &e){
        std::cout << e.what() << std::endl;
        throw;
    }
    catch(fs::filesystem_error &e){
        std::cout << e.what() << std::endl;
        throw IOError{e.what()};
    }
}
