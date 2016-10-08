#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

static PyObject* set_bonds(PyObject *self, PyObject *args)
{
    PyArrayObject *coord_py;
    PyArrayObject *cutoff_py;

    float cutfac;

    PyObject *off_py;

    float dist_n,dist_v[3];
    PyArray_Descr *float32 = PyArray_DescrFromType(NPY_FLOAT32);
    float effcut;

    if (!PyArg_ParseTuple(args, "OOfO", &coord_py, &cutoff_py, &cutfac, &off_py)) return NULL;
    if (!(PyArray_Check(coord_py)|PyArray_ISFLOAT(coord_py))) return NULL;
    if (!(PyArray_Check(cutoff_py)|PyArray_ISFLOAT(cutoff_py))) return NULL;
    if (!(PyList_Check(off_py)|(PyList_Size(off_py)==8))) return NULL;

    int nat = PyArray_DIM(cutoff_py,0);
    float *cutoff = PyArray_GETPTR1(cutoff_py,0);
    float *coord = PyArray_GETPTR2(coord_py,0,0);
    PyObject *bonds = PyList_New(8);
    for (Py_ssize_t i=0;i<8;i++){
        PyList_SetItem(bonds,i,PyList_New(0));
    }

    for (Py_ssize_t dir_i=0;dir_i<8;dir_i++) {
        PyObject *dir_py = PyList_GetItem(off_py,dir_i);
        PyObject *bond_off = PyList_GetItem(bonds,dir_i);

        for (Py_ssize_t off_i=0; off_i<PyList_Size(dir_py); off_i++) {
            PyObject *offset_py = PyList_GetItem(dir_py,off_i);
            float* off1 = PyArray_GETPTR1((PyArrayObject*) PyTuple_GetItem(offset_py,0),0);
            float* off2 = PyArray_GETPTR1((PyArrayObject*) PyTuple_GetItem(offset_py,1),0);

            for (int i=0;i<nat;i++) {
                if (cutoff[i] == 0) continue;

                for (int j=0;j<nat;j++) {
                    if ((j<i)&&(!memcmp(off1,off2,3*sizeof(float)))) continue;
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
                        dist_n=sqrtf(dist_n);
                        PyObject* distance = PyArray_Scalar(&dist_n,float32,NULL);
                        PyObject* bond = Py_BuildValue("iiOO",i,j,offset_py,distance);
                        PyList_Append(bond_off,bond);
                        Py_DECREF(bond);
                    }
                }
            }
        }
    }

    return bonds;
}

static PyObject* make_vol_gradient(PyObject *self, PyObject *args)
{
    PyArrayObject *volume_py;
    if (!PyArg_ParseTuple(args,"O",&volume_py)) return NULL;
    if (!PyArray_Check(volume_py)) return NULL;
    int x = PyArray_DIM(volume_py,0);
    int y = PyArray_DIM(volume_py,1);
    int z = PyArray_DIM(volume_py,2);
    npy_intp d[] = {3,x,y,z};
    int il,ih,jl,jh,kl,kh;
    PyArrayObject *gradient_py = (PyArrayObject*)PyArray_SimpleNew(4,d,PyArray_TYPE(volume_py));
    float *gradient = PyArray_GETPTR4(gradient_py,0,0,0,0);
    float *volume = PyArray_GETPTR3(volume_py,0,0,0);
    for (int i=0;i<x;i++) {
        if(i==0){
            il=x-1;ih=1;
        }else if(i==x-1){
            il=x-2;ih=0;
        }else{
            il=i-1;ih=i+1;
        }
        for (int j=0;j<y;j++) {
            if(j==0){
                jl=y-1;jh=1;
            }else if(j==y-1){
                jl=y-2;jh=0;
            }else{
                jl=j-1;jh=j+1;
            }
            for (int k=0;k<z;k++) {
                if(k==0){
                    kl=z-1;kh=1;
                }else if(k==z-1){
                    kl=z-2;kh=0;
                }else{
                    kl=k-1;kh=k+1;
                }
                gradient[0*x*y*z+i*y*z+j*z+k]=volume[il*y*z+j*z+k]-volume[ih*y*z+j*z+k];
                gradient[1*x*y*z+i*y*z+j*z+k]=volume[i*y*z+jl*z+k]-volume[i*y*z+jh*z+k];
                gradient[2*x*y*z+i*y*z+j*z+k]=volume[i*y*z+j*z+kl]-volume[i*y*z+j*z+kh];
            }
        }
    }
    return (PyObject*)gradient_py;
}

static PyMethodDef MolMethods[] = {
    {"setBondsC", set_bonds, METH_VARARGS,"Calculate bonds inside the cell"},
    {"makeVolGradient", make_vol_gradient, METH_VARARGS,"Calculate gradient of grid-data"},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >=3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mol_c",
    NULL,
    -1,
    MolMethods
};

PyObject * PyInit_mol_c(void)
#else
PyMODINIT_FUNC initmol_c(void)
#endif
{
    import_array();
#if PY_MAJOR_VERSION >=3
    PyObject *module = PyModule_Create(&moduledef);
    if (module == NULL) return NULL;
    return module;
#else
    Py_InitModule("mol_c", MolMethods);
#endif
}
