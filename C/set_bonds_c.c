#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

static PyObject* set_bonds(PyObject *self, PyObject *args)
{
    PyObject *coord_py;
    PyObject *cutoff_py;
    float **coord, *cutoff;
    float cutfac;
    PyObject *off1_py, *off2_py;
    float *off1, *off2;
    float dist_n,dist_v[3];
    float effcut;
    int nat;
    int status=0;

    if (!PyArg_ParseTuple(args, "OOf(OO)", &coord_py, &cutoff_py, &cutfac, &off1_py, &off2_py)) return NULL;
    /*TODO: check PyArrays for dimension and type*/
    if (!PyArray_Check(coord_py)) return NULL;
    if (!PyArray_Check(cutoff_py)) return NULL;
    if (!PyArray_Check(off1_py)) return NULL;
    if (!PyArray_Check(off2_py)) return NULL;

    nat = PyArray_DIMS((PyArrayObject*) cutoff_py)[0];
    cutoff = (float*)PyArray_GETPTR1((PyArrayObject*) cutoff_py,0);
    coord = (float**)PyArray_GETPTR2((PyArrayObject*) coord_py,0,0);
    off1 = (float*)PyArray_GETPTR1((PyArrayObject*) off1_py,0);
    off2 = (float*)PyArray_GETPTR1((PyArrayObject*) off2_py,0);

    for (int i=0;i<nat;i++) {
        if (cutoff[i] == 0) continue;
        for (int j=0;j<nat;j++) {
            printf("hallo!\n");
            if ((j<i)&&(memcmp(off1,off2,3*sizeof(float)))) continue;
            if (i==j) continue;
            if (cutoff[j]==0) continue;
            effcut=(cutoff[i]+cutoff[j])*cutfac;
            dist_v[0] = coord[i][0] + off1[0] - coord[j][0] - off2[0];
            if (dist_v[0]>effcut) continue;
            dist_v[1] = coord[i][1] + off1[1] - coord[j][1] - off2[1];
            if (dist_v[1]>effcut) continue;
            dist_v[2] = coord[i][2] + off1[2] - coord[j][2] - off2[2];
            if (dist_v[2]>effcut) continue;
            dist_n=dist_v[0]*dist_v[0]+dist_v[1]*dist_v[1]+dist_v[2]*dist_v[2];
            if((0.57<dist_n)&&(dist_n<effcut*effcut)){
                printf("bond of length %f between atoms %i and %i\n",dist_n,i,j);
            }else{
                printf("no bond between atoms %i and %i\n",i,j);
            }
        }
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
