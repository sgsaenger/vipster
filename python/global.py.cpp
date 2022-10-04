#include "global.py.h"
#include "atom.py.h"
#include "bond.py.h"
#include "data.py.h"
#include "configfile.py.h"
#include "fileio.py.h"
#include "parameters.py.h"
#include "plugin.py.h"
#include "presets.py.h"
#include "kpoints.py.h"
#include "molecule.py.h"
#include "periodictable.py.h"
#include "step.py.h"
#include "vec.py.h"

using namespace Vipster;

PYBIND11_MAKE_OPAQUE(std::map<std::string, std::string>)
PYBIND11_MAKE_OPAQUE(std::vector<std::string>)

void Vipster::Py::setupVipster(py::module &m, ConfigState &state, bool enableWrite)
{
    m.doc() = "Vipster\n"
              "=======\n\n"
              "A molecular modeling framework with periodic structures in mind.\n"
              "Use readFile() and writeFile() to handle files.\n"
              "Please inspect Molecule and Step as the main data "
              "containers for more information.";

    /*
     * Basic containers
     */

    py::bind_map<std::map<std::string,std::string>>(m, "__StrStrMap");
    py::bind_vector<std::vector<std::string>>(m, "__StrVector");
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
    Py::Data(m);
    // Read config state, init state-dependent API
    Py::Molecule(m, state);
    Py::FileIO(m, state, enableWrite);
    Py::Plugins(m, state);
    Py::Parameters(m);
    Py::Presets(m);
    // expose state
    Py::config(m, state);
}
