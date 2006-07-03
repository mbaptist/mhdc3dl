#include <Python.h>

#include "lss.h"

#include <string>

#include <iostream>

#include <cat.h>


#include <goops.h>

using namespace cat;

using namespace std;

//Launch the lss code
PyObject *
mhdc3dl_python_lss_run(PyObject *self, PyObject *args)
{
  double lambda_min,lambda_max;
  char * input_module_name;
  if (!PyArg_ParseTuple(args, "s", &input_module_name ))
    return NULL; 
	input input_obj(PyImport_ImportModule(input_module_name));
	lss lss_obj(input_obj);
	lss_obj.run(lambda_min,lambda_max);
  return Py_BuildValue("dd",lambda_min,lambda_max);
};

//Add Methods to function
PyMethodDef mhdc3dl_python_Methods[] = {
    {"lss_run", mhdc3dl_python_lss_run , METH_VARARGS,
     "Run large scale stability for a magnetoconvective system"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initmhdc3dl_python(void)
{
	(void) Py_InitModule("mhdc3dl_python", mhdc3dl_python_Methods);
};

