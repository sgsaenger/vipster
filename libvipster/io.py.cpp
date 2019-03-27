#include "pyvipster.h"
#include "io.h"

namespace Vipster::Py{
void PWInput(py::module&);
void LmpInput(py::module&);
void XYZ(py::module&);
void CPInput(py::module&);
void ORCA(py::module&);

void IO(py::module& m){
    auto io = m.def_submodule("IO");

    py::enum_<IOFmt>(io,"Fmt")
        .value("XYZ", IOFmt::XYZ)
        .value("PWI", IOFmt::PWI)
        .value("PWO", IOFmt::PWO)
        .value("LMP", IOFmt::LMP)
        .value("DMP", IOFmt::DMP)
        .value("CPI", IOFmt::CPI)
        .value("CUBE", IOFmt::CUBE)
        .value("XSF", IOFmt::XSF)
    ;

    m.def("readFile",[](std::string fn, IOFmt fmt)->py::object{
        auto data = readFile(fn,fmt);
        if(data.data.empty()){
            return py::make_tuple(data.mol, std::move(data.param), py::none());
        }else{
            py::list l{};
            for(auto& d: data.data){
                l.append(d.release());
            }
            return py::make_tuple(data.mol, std::move(data.param), l);
        }
    }, "filename"_a, "format"_a);

    /*
     * TODO: provide wrapper
     *
     * fall back to default-config/param
     * only index of state is relevant when not wrapping the GUI
     */
    py::class_<IO::State>(io, "State")
        .def_readwrite("index", &IO::State::index)
        .def_readwrite("atom_fmt", &IO::State::atom_fmt)
        .def_readwrite("cell_fmt", &IO::State::cell_fmt)
    ;

    m.def("writeFile", &writeFile, "filename"_a, "format"_a, "molecule"_a,
          "param"_a=nullptr, "config"_a=nullptr, "state"_a=IO::State{});

    py::class_<IO::BaseParam>(io, "BaseParam")
        .def_readwrite("name", &IO::BaseParam::name)
    ;

    py::class_<IO::BaseConfig>(io, "BaseConfig")
        .def_readwrite("name", &IO::BaseConfig::name)
    ;

    /*
     * Initialize plugins if needed
     */
    PWInput(io);
    LmpInput(io);
    XYZ(io);
    CPInput(io);
    ORCA(io);
}
}
