#include <Python.h>

#include "lss.h"

#include <string>

#include <iostream>

#include <cat.h>


#include <goops.h>

using namespace cat;

using namespace std;



//Launch the lss code
static PyObject *
mhdc3dl_python_lss_run(PyObject *self, PyObject *args)
{

  double lambda_minimal,lambda_maximal;
  char * input_module_name;
  if (!PyArg_ParseTuple(args, "s", &input_module_name ))
    return NULL;
  
  input input_obj(PyImport_ImportModule(input_module_name));
 
  lss lss_obj(input_obj);
  lss_obj.run(lambda_minimal,lambda_maximal);

  //lambda_minimal=0;
  //lambda_maximal=1;

  return Py_BuildValue("dd",lambda_minimal,lambda_maximal);
};

//Add Methods to function
static PyMethodDef mhdc3dl_python_Methods[] = {
    {"lss_run", mhdc3dl_python_lss_run , METH_VARARGS,
     "Run large scale stability for a magnetoconvective system"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initmhdc3dl_python(void)
{
    (void) Py_InitModule("mhdc3dl_python", mhdc3dl_python_Methods);
};

