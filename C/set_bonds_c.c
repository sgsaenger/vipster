#include <python/Python.h>
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

static PyObject* set_bonds(PyObject *self, PyObject *args)
{
    PyObject *coord;
    PyObject *cutoff;
    float cutfac;
    PyObject *off1, *off2;
    int status;

    if(!PyArg_ParseTuple(args, "OOf(OO)", &coord, &cutoff, &cutfac, &off1, &off2)) return NULL;
    //if(!PyArg_ParseTuple(args,"O",&coord)) return NULL;
    status = PyArray_Check(coord);
    for(int i=0;i<PyArray_NDIM(coord);i++){
        printf("%i\n",PyArray_DIMS(coord)[i]);
    }
    status = PyArray_Check(cutoff);
    for(int i=0;i<PyArray_NDIM(cutoff);i++){
        printf("%i\n",PyArray_DIMS(cutoff)[i]);
    }
    status = PyArray_Check(off1);
    for(int i=0;i<PyArray_NDIM(off1);i++){
        printf("%i\n",PyArray_DIMS(off1)[i]);
    }
    status = PyArray_Check(off2);
    for(int i=0;i<PyArray_NDIM(off2);i++){
        printf("%i\n",PyArray_DIMS(off2)[i]);
    }

    return Py_BuildValue("i",status);
}

static PyMethodDef BondsMethods[] = {
    {"set_bonds_c", set_bonds, METH_VARARGS,"Execute a shell command."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initbonds(void)
{
    import_array();
    Py_InitModule("bonds", BondsMethods);
}
