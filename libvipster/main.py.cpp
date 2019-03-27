#include "pyvipster.h"

#include <sstream>
#include "molecule.h"
#include "io.h"

using namespace Vipster;

namespace Vipster{
namespace Py{
void Vec(py::module&);
void Atom(py::module&);
void Bond(py::module&);
void PSE(py::module&);
void Step(py::module&);
void Selection(py::module&);
void KPoints(py::module&);
void Molecule(py::module&);
void IO(py::module&);
void Data(py::module&);
}
}

PYBIND11_MODULE(vipster, m) {
    m.doc() = "Vipster\n"
              "=======\n\n"
              "A molecular modeling framework with periodic structures in mind.\n"
              "Use readFile() and writeFile() to handle files.\n"
              "Please inspect Molecule and Step as the main data"
              "containers for more information.";

    /*
     * Basic containers
     */

    py::bind_map<std::map<std::string,std::string>>(m, "__StrStrMap__");
    py::bind_vector<std::vector<std::string>>(m, "__StrVector__");
    bind_array<ColVec>(m, "ColVec");

    /*
     * Initialize library
     */

    Py::Vec(m);
    Py::Atom(m);
    Py::Bond(m);
    Py::PSE(m);
    Py::Step(m);
    Py::KPoints(m);
    Py::Molecule(m);
    Py::IO(m);
    Py::Data(m);
}
