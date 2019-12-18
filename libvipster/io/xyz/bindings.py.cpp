#include "pyvipster.h"
#include "plugin.h"

namespace Vipster::Py{
void XYZ(py::module& m){
//    auto c = py::class_<IO::XYZPreset, IO::BasePreset>(m, "XYZPreset");

//    py::enum_<IO::XYZPreset::Data>(c, "Data")
//        .value("None", IO::XYZPreset::Data::None)
//        .value("Charge", IO::XYZPreset::Data::Charge)
//        .value("Forces", IO::XYZPreset::Data::Forces)
//    ;

//    py::enum_<IO::XYZPreset::Mode>(c, "Mode")
//        .value("Step", IO::XYZPreset::Mode::Step)
//        .value("Trajec", IO::XYZPreset::Mode::Trajec)
//        .value("Cell", IO::XYZPreset::Mode::Cell)
//    ;

//    c.def(py::init<IO::XYZPreset::Mode, IO::XYZPreset::Data>(),
//          "filemode"_a=IO::XYZPreset::Mode::Step,
//          "atomdata"_a=IO::XYZPreset::Data::None)
//     .def_readwrite("filemode", &IO::XYZPreset::filemode)
//     .def_readwrite("atomdata", &IO::XYZPreset::atomdata)
//    ;
}
}
