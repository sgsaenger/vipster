#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>

static PyObject* set_bonds(PyObject *self, PyObject *args)
{
    PyObject *coord_py;
    PyObject *cutoff_py;
    float *coord, *cutoff;
    int nat;

    float cutfac;

    PyObject *off_py;

    float dist_n,dist_v[3];
    float effcut;

    PyObject *bonds;

    if (!PyArg_ParseTuple(args, "OOfO", &coord_py, &cutoff_py, &cutfac, &off_py)) return NULL;
    if (!(PyArray_Check(coord_py)|PyArray_ISFLOAT(coord_py))) return NULL;
    if (!(PyArray_Check(cutoff_py)|PyArray_ISFLOAT(cutoff_py))) return NULL;
    if (!(PyList_Check(off_py)|(PyList_Size(off_py)==8))) return NULL;

    nat = PyArray_DIMS((PyArrayObject*) cutoff_py)[0];
    cutoff = (float*)PyArray_GETPTR1((PyArrayObject*) cutoff_py,0);
    coord = (float*)PyArray_GETPTR2((PyArrayObject*) coord_py,0,0);
    bonds = PyList_New(8);

    for (Py_ssize_t dir_i=0;dir_i<3;dir_i++) {
        PyObject *dir_py = PyList_GetItem(off_py,dir_i);
        PyObject *bond_off = PyList_New(0);

        for (Py_ssize_t off_i=0; off_i<PyList_Size(dir_py); off_i++) {
            PyObject *offset_py = PyList_GetItem(dir_py,off_i);
            float* off1 = PyArray_GETPTR1((PyArrayObject*) PyTuple_GetItem(offset_py,0),0);
            float* off2 = PyArray_GETPTR1((PyArrayObject*) PyTuple_GetItem(offset_py,1),0);

            for (int i=0;i<nat;i++) {
                if (cutoff[i] == 0) continue;

                for (int j=0;j<nat;j++) {
                    if ((j<i)&&(memcmp(off1,off2,3*sizeof(float)))) continue;
                    if (i==j) continue;
                    if (cutoff[j]==0) continue;
                    effcut=(cutoff[i]+cutoff[j])*cutfac;
                    dist_v[0] = coord[i*3+0] - coord[j*3+0];
                    dist_v[0] = coord[i*3+0] + off1[0] - coord[j*3+0] - off2[0];
                    if (dist_v[0]>effcut) continue;
                    dist_v[1] = coord[i*3+1] + off1[1] - coord[j*3+1] - off2[1];
                    if (dist_v[1]>effcut) continue;
                    dist_v[2] = coord[i*3+2] + off1[2] - coord[j*3+2] - off2[2];
                    if (dist_v[2]>effcut) continue;
                    dist_n=dist_v[0]*dist_v[0]+dist_v[1]*dist_v[1]+dist_v[2]*dist_v[2];
                    if((0.57<dist_n)&&(dist_n<effcut*effcut)){
                        PyObject* bond = Py_BuildValue("iiOf",i,j,offset_py,sqrt(dist_n));
                        PyList_Append(bond_off,bond);
                        Py_DECREF(bond);
                    }
                }
            }
        }
        PyList_SetItem(bonds,dir_i,bond_off);
    }

    return bonds;
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
