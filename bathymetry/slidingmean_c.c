/*  Fast C implementation of sliding mean
    Usage: df = slidingmean_c(t,d,w)
    Where: t,d are 1D numpy arrays, type float64, length N
           w  is a 1D numpy array,  type float64, length 1
    Michael Dunphy Nov 2016
*/

#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>

void smean(double *t, double *d, double w, double *sa, size_t n) {
    // input vectors t and d of length n, output sa also length n
    size_t i,j,is=0;
    for(i=1;i<n;i++)
    {
        if (t[i]-t[0] < w) { /* in the first w */
            sa[i]=0.0;
            for(j=0;j<i;j++)
                sa[i] += 0.5*(d[j]+d[j+1])*(t[j+1]-t[j]);
            sa[i] /= t[i]-t[0];
        }
        else {
            while (t[i]-t[is] > w) is++;
            /* if (is==0) printf("Problem is=%d, should be >0\n",is); */
            sa[i] = 0.5*(d[is-1]+d[is])*(t[is]-t[is-1]);
            sa[i] *= (t[is]-t[i]+w) / (t[is]-t[is-1]);
            for (j=is;j<i;j++)
                    sa[i] += 0.5*(d[j]+d[j+1])*(t[j+1]-t[j]);
            sa[i] /= w;
        }
    }
}

// Main slidingmean_c function
static PyObject* slidingmean_c(PyObject* self, PyObject* args)
{
    PyArrayObject *pt, *pd, *pw, *pr=NULL;
    double        *t, *d, *w, *r, width;
    npy_uint      N;

    // Parse numpy ndarray arguments
    if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &pt,
        &PyArray_Type, &pd, &PyArray_Type, &pw ))
        return NULL;

    npy_uint siz1 = PyArray_ITEMSIZE(pt);
    npy_uint siz2 = PyArray_ITEMSIZE(pd);
    if ((siz1 != sizeof(double)) || (siz2 != sizeof(double)) )
        printf("Expecting 8 byte double data, given %d and %d byte data\n",siz1,siz2);

    npy_uint chk1 = PyArray_CHKFLAGS(pt, NPY_ARRAY_C_CONTIGUOUS);
    npy_uint chk2 = PyArray_CHKFLAGS(pd, NPY_ARRAY_C_CONTIGUOUS);
    if (!chk1 || !chk2) printf("Expecting C ordered arrays\n");

    // Create empty output array
    npy_intp *dims = PyArray_DIMS(pd);
    pr = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);

    // Set up pointers to input
    t = (double *) PyArray_DATA(pt);
    d = (double *) PyArray_DATA(pd);
    w = (double *) PyArray_DATA(pw); width = w[0];
    r = (double *) PyArray_DATA(pr);
    N = PyArray_DIM(pt,0);

    // Call the sliding mean function
    smean(t,d,width,r,N);

    return PyArray_Return(pr);
}


// Define methods (same for Python 2 and 3)
static PyMethodDef slidingmean_c_methods[] =
{
    {"slidingmean_c", slidingmean_c, METH_VARARGS, "Fast C implementation of sliding mean." }, {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION < 3
// Module initialization for Python 2
PyMODINIT_FUNC initslidingmean_c(void)
{
    (void) Py_InitModule("slidingmean_c", slidingmean_c_methods);
    import_array();
}

#else
// Define module for Python 3
static struct PyModuleDef slidingmean_c_definition = {
    PyModuleDef_HEAD_INIT,
    "slidingmean_c",
    "C implementation of sliding mean.",
    -1,
    slidingmean_c_methods
};

// Module initialization for Python 3
PyMODINIT_FUNC PyInit_slidingmean_c(void)
{
    Py_Initialize();
    PyObject *m = PyModule_Create(&slidingmean_c_definition);
    import_array();
    return m;
}
#endif

