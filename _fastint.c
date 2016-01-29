//
//  _fastint.c
//  PythonProteins
//
//  Created by Venkatesh Sivaraman on 12/19/14.
//  Copyright (c) 2014 Venkatesh Sivaraman. All rights reserved.
//

#include <Python.h>

/* Docstrings */
static char fastint_docstring[] = "Quickly convert a string to an integer.";

/* Available functions */
static PyObject *fastint_int(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
	{"fastint", fastint_int, METH_VARARGS, fastint_docstring},
	{NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC initfastint(void)
{
	PyObject *m = Py_InitModule("fastint", module_methods);
	if (m == NULL)
		return;
}

static PyObject *fastint_int(PyObject *self, PyObject *args) {
	char *s; unsigned r = 0;
	if (!PyArg_ParseTuple(args, "s", &s)) return NULL;
	for (r = 0; *s; r = r * 10 + *s++ - '0');
	return Py_BuildValue("i", r);
}