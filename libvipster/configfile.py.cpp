#include "pyvipster.h"
#include "configfile.h"
#include "io/plugin.h"

namespace Vipster::Py{
void config(py::module& m, const ConfigState& state){
    m.attr("PeriodicTable") = std::get<0>(state);
//    m.attr("Plugins") = std::move(std::get<2>(state));
//    m.attr("Parameters") = std::move(std::get<3>(state));
//    m.attr("IOPresets") = std::move(std::get<4>(state));
    m.attr("Plugins") = py::cast(std::get<2>(state), py::return_value_policy::reference);
    m.attr("Parameters") = py::cast(std::get<3>(state), py::return_value_policy::reference);
    m.attr("IOPresets") = py::cast(std::get<4>(state), py::return_value_policy::reference);
//    m.attr("IOPresets") = std::get<4>(state);

//    m.def("__readConfig", [&](){
//        auto [t, s, p, c] = readConfig();
//        return std::move(p.at(IOFmt::PWI));
//        auto object = py::make_tuple<py::return_value_policy::move>(t, s,
//                std::move(p));//, std::move(c));
//        m.attr("PeriodicTable") = std::get<0>(tmp);
//        py::print("hallo");
//        m.attr("Settings").cast<Settings>();// = Settings{};//std::get<1>(tmp);
//        py::print("was?");
//        m.attr("Parameters") = std::move(std::get<2>(tmp));
//        m.attr("IOPresets") = std::move(std::get<3>(tmp));
//    });
//    }, py::return_value_policy::move);

//    m.def("__saveConfig", [&](){
//        auto tmp = readConfig();

        /* try to use module variables
         * if they don't match the type, fall back to defaults
         */
//        PeriodicTable tmpTable = std::get<0>(tmp);
//        try{ tmpTable = m.attr("PeriodicTable").cast<PeriodicTable>(); } catch(py::cast_error){}
//        Settings tmpSet = std::get<1>(tmp);
//        try{ tmpSet = m.attr("Settings").cast<Settings>(); } catch(py::cast_error){}
//        IO::Parameters tmpParam = std::move(std::get<2>(tmp));
//        try{ tmpParam = m.attr("Parameters").cast<IO::Parameters>(); } catch(py::cast_error){}
//        IO::Configs tmpPreset = std::move(std::get<3>(tmp));
//        try{ tmpPreset = m.attr("IOPresets").cast<IO::Configs>(); } catch(py::cast_error){}

//        ConfigState newConf{ tmpTable, tmpSet, tmpParam, tmpPreset};
//        try{
//            saveConfig(newConf);
//        }catch(IO::Error&e){
//            return e.what();
//        }
//        return "";
//    });
}
}
