#include "pyvipster.h"
#include "data.h"

namespace Vipster::Py{
void Data(py::module& m){
   auto b = py::class_<BaseData>(m, "BaseData")
       .def_readwrite("name", &BaseData::name)
   ;
   auto f2 = py::class_<DataGrid2D_f, BaseData>(m, "DataGrid2D_f")
       .def_readwrite("cell", &DataGrid2D_f::cell)
       .def_readwrite("origin", &DataGrid2D_f::origin)
       .def_readonly("extent", &DataGrid2D_f::extent)
       .def_readonly("size", &DataGrid2D_f::size)
   ;
   bind_array<DataGrid2D_f::Extent>(f2, "Extent");
   auto f3 = py::class_<DataGrid3D_f, BaseData>(m, "DataGrid3D_f")
       .def_readwrite("cell", &DataGrid3D_f::cell)
       .def_readwrite("origin", &DataGrid3D_f::origin)
       .def_readonly("extent", &DataGrid3D_f::extent)
       .def_readonly("size", &DataGrid3D_f::size)
   ;
   bind_array<DataGrid3D_f::Extent>(f3, "Extent");
}
}
