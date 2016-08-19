#include "xyz.h"
#include <fstream>

using namespace Vipster;

std::tuple<Molecule,optional<Param>> Vipster::IO::XYZ_parser(std::string fn, std::ifstream &file)
{
    m = Molecule(fn,0);
    return std::tuple<Molecule,optional<Param>>(m,optional<Param>());
}

void    Vipster::IO::XYZ_writer(const Molecule &m, std::ofstream &file,Param p)
{

}
