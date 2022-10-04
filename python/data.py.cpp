#include "data.py.h"
#include "vipster/data.h"

template <typename T>
py::class_<T, Vipster::BaseData> bind_datagrid(py::module &m, std::string const &name){
   auto tmp = py::class_<T, Vipster::BaseData>{m, name.c_str()}
       .def(py::init([](py::args args){
           py::tuple obj;
           if((args.size() == 1) &&
              ((py::isinstance<py::tuple>(args[0]) ||
               (py::isinstance<py::list>(args[0]))))){
                obj = py::cast<py::tuple>(args[0]);
           }else{
                    obj = args;
           }
           if(obj.size() != T::Dim){
               throw py::type_error();
           }
           typename T::Extent extent;
           for(size_t i=0; i<T::Dim; ++i){
               extent[i] = py::cast<size_t>(obj[i]);
           }
           return T{extent};
       }))
       .def_readwrite("cell", &T::cell)
       .def_readwrite("origin", &T::origin)
       .def_readonly("extent", &T::extent)
       .def_readonly("shape", &T::extent)
       .def_readonly("size", &T::size)
       .def("__getitem__", [](T& d, py::tuple& t){
            if(t.size() != d.Dim){
                throw py::type_error();
            }
            typename T::Extent idx;
            for(size_t i=0; i<T::Dim; ++i){
                auto j = py::cast<long>(t[i]);
                if(j < 0){
                    j = j + static_cast<long>(d.extent[i]);
                }
                if((j < 0) || j >= static_cast<long>(d.extent[i])){
                    throw py::index_error();
                }
                idx[i] = static_cast<size_t>(j);
            }
            return d(idx);
       })
       .def("__setitem__", [](T& d, py::tuple& t, typename T::value_type val){
            if(t.size() != d.Dim){
                throw py::type_error();
            }
            typename T::Extent idx;
            for(size_t i=0; i<T::Dim; ++i){
                auto j = py::cast<long>(t[i]);
                if(j < 0){
                    j = j + static_cast<long>(d.extent[i]);
                }
                if((j < 0) || j >= static_cast<long>(d.extent[i])){
                    throw py::index_error();
                }
                idx[i] = static_cast<size_t>(j);
            }
            d(idx) = val;
       })
       .def("__len__", [](const T &d){return d.extent[0];})
   ;
   return tmp;
}

void Vipster::Py::Data(py::module& m){
   auto b = py::class_<BaseData>(m, "__BaseData")
       .def_readwrite("name", &BaseData::name)
   ;

   bind_datagrid<DataGrid2D_f>(m, "DataGrid2D_f");
   bind_array<typename DataGrid2D_f::Extent>(b, "Extent2D");
   bind_datagrid<DataGrid3D_f>(m, "DataGrid3D_f");
   bind_array<typename DataGrid3D_f::Extent>(b, "Extent3D");
   bind_datagrid<DataGrid3D_v>(m, "DataGrid3D_v");
}
