#include "pyvipster.h"

#include <sstream>
#include "molecule.h"
#include "fileio.h"
#include "configfile.h"

using namespace Vipster;

namespace Vipster::Py{
static ConfigState state;
void Vec(py::module&);
void Atom(py::module&);
void Bond(py::module&);
void Table(py::module&);
void Step(py::module&);
void KPoints(py::module&);
void Molecule(py::module&);
void Data(py::module&);
void IO(py::module&, const ConfigState&);
void config(py::module&, const ConfigState&);
}

PYBIND11_MODULE(_vipster, m) {
    m.doc() = "Vipster\n"
              "=======\n\n"
              "A molecular modeling framework with periodic structures in mind.\n"
              "Use readFile() and writeFile() to handle files.\n"
              "Please inspect Molecule and Step as the main data "
              "containers for more information.";

    /*
     * Basic containers
     */

    py::bind_map<std::map<std::string,std::string>>(m, "__StrStrMap__");
    py::bind_vector<std::vector<std::string>>(m, "__StrVector__");
    bind_array<ColVec>(m, "ColVec");
    bind_array<DiffVec>(m, "DiffVec");

    /*
     * Initialize library
     */


    Py::Vec(m);
    Py::Atom(m);
    Py::Bond(m);
    Py::Table(m);
    Py::Step(m);
    Py::KPoints(m);
    Py::Molecule(m);
    Py::Data(m);
    // I/O
    Py::state = readConfig();
    Py::IO(m, Py::state);
    Py::config(m, Py::state);
}
