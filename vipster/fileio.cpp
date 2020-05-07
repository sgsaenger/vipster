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

const IO::Plugin *Vipster::guessFmt(std::string fn, const IO::Plugins &p)
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
        throw IO::Error{"Could not deduce format of file \""+fn+
                        "\"\nPlease specify format explicitely", false};
    }
    // read file
    return readFile(fn, plugin);
}

// read with explicit format
IO::Data Vipster::readFile(const std::string &fn, const IO::Plugin *plug)
{
    // set locale to C to get consistent parsing
    std::string userLocale = setlocale(0, nullptr);
    setlocale(LC_ALL, "C");
    // check if file can be read
    std::ifstream file{fn};
    if(!file){
        setlocale(LC_ALL, userLocale.c_str());
        throw IO::Error("Could not open \""+fn+'"');
    }
    // try to parse
    if(!plug->parser){
        setlocale(LC_ALL, userLocale.c_str());
        throw IO::Error("Format is not readable");
    }
    auto tmp = plug->parser(fn, file);
    if(!tmp.mol.getNstep()){
        setlocale(LC_ALL, userLocale.c_str());
        throw IO::Error("No Molecule could be parsed");
    }
    // return if successful
    setlocale(LC_ALL, userLocale.c_str());
    return tmp;
}

bool  Vipster::writeFile(const std::string &fn,
                         const IO::Plugin *plug,
                         const Molecule &m,
                         std::optional<size_t> idx,
                         const std::optional<IO::Parameter>& p,
                         const std::optional<IO::Preset>& c)
{
    if(!idx){
        idx = m.getNstep()-1;
    }
    try{
        if(!plug->writer){
            throw IO::Error{"Read-only format"};
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
                    throw IO::Error{"Could not open "+fn};
                }
            }
            res = plug->writer(m, file, p, c, *idx);
        }
        if(use_temp){
            if(fs::exists(fn)){
                fs::remove(fn); // overwrite_existing not working on windows as of mingw-w64 gcc9
            }
            fs::copy_file(filename, fn, fs::copy_options::overwrite_existing);
        }
        return res;
    }
    catch(IO::Error &e){
        std::cout << e.what() << std::endl;
        throw;
    }
    catch(fs::filesystem_error &e){
        std::cout << e.what() << std::endl;
        throw IO::Error{e.what()};
    }
}
