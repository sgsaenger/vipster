#include "pyvipster.h"
#include "io/xyz/plugin.h"

namespace Vipster::Py{
void XYZ(py::module& m){
    auto c = py::class_<IO::XYZConfig>(m, "XYZConfig")
        .def_readwrite("filemode", &IO::XYZConfig::filemode)
        .def_readwrite("atomdata", &IO::XYZConfig::atomdata)
    ;

    py::enum_<IO::XYZConfig::Data>(c, "Data")
        .value("None", IO::XYZConfig::Data::None)
        .value("Charge", IO::XYZConfig::Data::Charge)
        .value("Forces", IO::XYZConfig::Data::Forces)
    ;

    py::enum_<IO::XYZConfig::Mode>(c, "Mode")
        .value("Step", IO::XYZConfig::Mode::Step)
        .value("Trajec", IO::XYZConfig::Mode::Trajec)
        .value("Cell", IO::XYZConfig::Mode::Cell)
    ;
}
}
