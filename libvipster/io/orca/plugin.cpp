#include "plugin.h"
#include "../util.h"

using namespace Vipster;

IO::Data OrcaParser(const std::string& name, std::ifstream &file){
    IO::Data d{};
    d.param = std::make_unique<IO::OrcaParam>(name);
    auto& p = static_cast<IO::OrcaParam&>(*d.param.get());
    Molecule &m = d.mol;
    m.setName(name);
    auto s = &m.newStep();

    bool foundCoords{false};
    AtomFmt at_fmt{AtomFmt::Angstrom};
    std::string line, c_fmt;
    int charge;
    unsigned int multiplicity;
    while(std::getline(file, line)){
        std::stringstream ls{line};
        std::string test;
        // prepare parsers for specific formats
        const auto xyzParser = [](Step& s, std::ifstream &file, std::string term){
            std::string line, test;
            std::stringstream ls;
            Vec coord;
            while(std::getline(file, line)){
                ls = std::stringstream{line};
                ls >> test;
                if(test == term){
                    return;
                }else{
                    ls >> coord[0] >> coord[1] >> coord[2];
                    s.newAtom(test, coord);
                }
            }
            throw IO::Error("ORCA: unterminated coordinates");
        };
        const auto intParser = [](Step& s, std::ifstream &file, std::string term){
            std::string line, test;
            std::stringstream ls;
            SizeVec ids;
            Vec values;
            while(std::getline(file, line)){
                ls = std::stringstream{line};
                ls >> test;
                if(test == term){
                    return;
                }else{
                    ls >> ids[0] >> ids[1] >> ids[2];
                    ls >> values[0] >> values[1] >> values[2];
                    for(const auto& id: ids){
                        if(id > s.getNat()){
                            throw IO::Error("ORCA: internal position refers"
                                            "to not-yet parsed atom: "+line);
                        }
                    }
                    IO::intToCart(s, test, ids, values);
                }
            }
        };
        const auto gzmtParser = [](Step& s, std::ifstream &file, std::string term){
            std::string line, test;
            std::stringstream ls;
            SizeVec ids;
            Vec values;
            while(std::getline(file, line)){
                ls = std::stringstream{line};
                ls >> test;
                if(test == term){
                    return;
                }else{
                    for(size_t i=0; i<3; ++i){
                        ls >> ids[i] >> values[i];
                        if(!ls.good()){
                            ids[i] = 0; values[i] = 0;
                        }
                        if(ids[i] > s.getNat()){
                            throw IO::Error("ORCA: internal position refers"
                                            "to not-yet parsed atom: "+line);
                        }
                    }
                    IO::intToCart(s, test, ids, values);
                }
            }
        };
        // parse file
        ls >> test;
        if(test == "*"){
            /* *: conventional coordinate input
             *
             * prefixed by format, total charge and multiplicity
             *
             */
            if(foundCoords){
                throw IO::Error("ORCA: multiple coordinate definitions");
            }else{
                foundCoords = true;
            }
            ls >> c_fmt >> charge >> multiplicity;
            std::transform(c_fmt.begin(), c_fmt.end(), c_fmt.begin(), ::tolower);
            if(c_fmt == "xyz"){
                xyzParser(*s, file, "*");
            }else if(c_fmt.substr(0, 3) == "int"){
                intParser(*s, file, "*");
            }else if(c_fmt == "gzmt"){
                gzmtParser(*s, file, "*");
            }else{
                throw IO::Error("ORCA: unsupported format descriptor: "+c_fmt);
            }
        }else{
            std::transform(test.begin(), test.end(), test.begin(), ::tolower);
            if(test == "%coords"){
                /* %coords: block-style coordinate input
                 *
                 */
                if(foundCoords){
                    throw IO::Error("ORCA: multiple coordinate definitions");
                }else{
                    foundCoords = true;
                }
                while(true){
                    ls >> test;
                    if(ls.good()){
                        std::transform(test.begin(), test.end(), test.begin(), ::tolower);
                        if(test == "ctyp"){
                            ls >> c_fmt;
                        }else if(test == "charge"){
                            ls >> charge;
                        }else if(test == "mult"){
                            ls >> multiplicity;
                        }else if(test == "units"){
                            ls >> test;
                            std::transform(test.begin(), test.end(), test.begin(), ::tolower);
                            if(test == "bohrs") at_fmt = AtomFmt::Bohr;
                        }else if(test == "coords"){
                            if(c_fmt == "xyz"){
                                xyzParser(*s, file, "end");
                            }else if(c_fmt == "internal"){
                                intParser(*s, file, "*");
                            }else{
                                throw IO::Error("ORCA: unsupported format descriptor in %coords-block: "+c_fmt);
                            }
                        }else if(test == "end"){
                            break;
                        }else if(test[0] == '#'){
                            continue;
                        }else{
                            throw IO::Error("ORCA: unknown line in %coords-block: "+line);
                        }
                    }else{
                        if(!std::getline(file, line)){
                            throw IO::Error("ORCA: reached end-of-file during parsing of %coords-block");
                        }
                        ls = std::stringstream{line};
                    }
                }
            }else if(!foundCoords){
                /* Header
                 *
                 * contains keword lines (starts with '!'),
                 * input blocks (starting with '%', ending with "end"),
                 * memory infos ("%maxcore N"),
                 * interspersed with comments ('#' until '#' or '\n')
                 *
                 * for now, parsed and saved as monolithic block
                 *
                 */
                if(test[0] == '!'){
                    ls >> test;
                    std::transform(test.begin(), test.end(), test.begin(), ::tolower);
                    if(test == "bohrs") at_fmt = AtomFmt::Bohr;
                }
                p.header.push_back(line);
            }//else{
                // TODO
                /* $new_job: trigger new job
                 *
                 * below this line, all previous sections may appear anew
                 *
                 * will create a new Vipster::Step
                 *
                 * header will be ignored
                 */
            //}
        }
    }
    s->modScale(at_fmt);

    return d;
}

bool OrcaWriter(const Molecule& m, std::ofstream &file,
                const IO::BaseParam *const p,
                const IO::BaseConfig *const,
                IO::State state)
{
    const auto *pp = dynamic_cast<const IO::OrcaParam*>(p);
    if(!pp) throw IO::Error("OrcaWriter needs ORCA parameter set");

    auto af = AtomFmt::Angstrom; // TODO: deduce from parameter set

    const auto& s = m.getStep(state.index).asFmt(af);

    for(const auto& line: pp->header){
        file << line << '\n';
    }
    // TODO: multiplicity + global charge
    file << "* xyz 0 1\n";
    for(const auto& at: s){
        file << std::left << std::setw(3) << at.name << " "
             << std::right << std::setw(10) << at.coord[0] << " "
             << std::right << std::setw(10) << at.coord[1] << " "
             << std::right << std::setw(10) << at.coord[2] << '\n';
    }
    file << "*\n";

    return true;
}

static std::unique_ptr<IO::BaseParam> makeParam(const std::string& name)
{
    return std::make_unique<IO::OrcaParam>(name);
}

const IO::Plugin IO::OrcaInput =
{
    "ORCA Input File",
    "inp",
    "orca",
    IO::Plugin::Param,
    &OrcaParser,
    &OrcaWriter,
    &makeParam,
    nullptr
};
