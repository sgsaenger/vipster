#include "global.py.h"

#include <sstream>

using namespace Vipster;

namespace Vipster::Py{
static ConfigState state;
}

PYBIND11_MODULE(vipster, m) {
    Py::state = readConfig();
    Py::setupVipster(m, Py::state, true);
}
