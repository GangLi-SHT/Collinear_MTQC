#include <Python.h>
#include <assert.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "mspglib.h"

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * py_get_dataset(PyObject *self, PyObject *args);
static PyObject * py_get_error_message(PyObject *self, PyObject *args);

struct module_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
  struct module_state *st = GETSTATE(m);
  PyErr_SetString(st->error, "something bad happened");
  return NULL;
}

static PyMethodDef _mspglib_methods[] = {
  {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
  {"dataset", py_get_dataset, METH_VARARGS, "Dataset for magnetic crystal symmetry"},
  {"error_message", py_get_error_message, METH_VARARGS, "Error message"},

  {NULL,NULL,0,NULL}
};

/*initialize*/

#if PY_MAJOR_VERSION >= 3

static int _mspglib_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int _mspglib_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_mspglib",
  NULL,
  sizeof(struct module_state),
  _mspglib_methods,
  NULL,
  _mspglib_traverse,
  _mspglib_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__mspglib(void)

#else
#define INITERROR return

  void
  init_mspglib(void)
#endif
{
  struct module_state *st;
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  PyObject *module = Py_InitModule("_mspglib", _mspglib_methods);
#endif

  if (module == NULL)
    INITERROR;

  st = GETSTATE(module);

  st->error = PyErr_NewException("_mspglib.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}

/*
static PyObject * py_get_version(PyObject *self, PyObject *args)
{
  PyObject *array;
  int i;
  int version[3];

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  version[0] = spg_get_major_version();
  version[1] = spg_get_minor_version();
  version[2] = spg_get_micro_version();

  array = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(array, i, PyLong_FromLong((long)version[i]));
  }

  return array;
}
*/

static PyObject * py_get_dataset(PyObject *self, PyObject *args)
{ 
  /*define input variable*/
  double symprec;
  PyArrayObject* py_lattice;
  PyArrayObject* py_positions;
  PyArrayObject* py_magmom;
  PyArrayObject* py_atom_types;

  PyObject *array, *vec, *mat, *rot, *trans, *time_reversal;
  PyObject *primitive_lattice, *mapping_to_primitive, *primitive_types;
  PyObject *positions, *magmom;


  int i, j, k, n;
  double (*lat)[3];
  double (*pos)[3];
  double (*mag)[3];
  int num_atom, len_list;
  int* typat;
  MspglibDataset *dataset;


  if (!PyArg_ParseTuple(args, "OOOOd",
                        &py_lattice,
                        &py_positions,
                        &py_magmom,
                        &py_atom_types,
                        &symprec)) {
    return NULL;
  }
  lat = (double(*)[3])PyArray_DATA(py_lattice);
  pos = (double(*)[3])PyArray_DATA(py_positions);
  mag = (double(*)[3])PyArray_DATA(py_magmom);
  num_atom = PyArray_DIMS(py_positions)[0];
  typat = (int*)PyArray_DATA(py_atom_types);


  if ((dataset = mspg_get_magnetic_dataset(lat,
                                           pos,
                                           mag,
                                           typat,
                                           num_atom,
                                           symprec)) == NULL) {
    Py_RETURN_NONE;
  }

  len_list = 17;
  array = PyList_New(len_list);
  n = 0;
  /* Space group number, international symbol, hall symbol */
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->spacegroup_number));
  n++;
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->hall_number));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->international_symbol));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->hall_symbol));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->choice));
  n++;
  
  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }
  PyList_SetItem(array, n, mat);
  n++;

  /* Origin shift */
  vec = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(vec, i, PyFloat_FromDouble(dataset->origin_shift[i]));
  }
  PyList_SetItem(array, n, vec);
  n++;
  
  /* Rotation matrices */
  rot = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    mat = PyList_New(3);
    for (j = 0; j < 3; j++) {
      vec = PyList_New(3);
      for (k = 0; k < 3; k++) {
        PyList_SetItem(vec, k, PyLong_FromLong((long) dataset->rotations[i][j][k]));
      }
      PyList_SetItem(mat, j, vec);
    }
    PyList_SetItem(rot, i, mat);
  }
  PyList_SetItem(array, n, rot);
  n++;


  /* Translation vectors */
  trans = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->translations[i][j]));
    }
    PyList_SetItem(trans, i, vec);
  }
  PyList_SetItem(array, n, trans);
  n++;

  /* Time reversal */
  time_reversal = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    PyList_SetItem(time_reversal, i, PyLong_FromLong((long) dataset->time_reversal[i])); 
  }
  PyList_SetItem(array, n, time_reversal);
  n++;



  primitive_lattice = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->primitive_lattice[i][j]));
    }
    PyList_SetItem(primitive_lattice, i, vec);
  }
  PyList_SetItem(array, n, primitive_lattice);
  n++;
  
  /* mapping table */
  mapping_to_primitive = PyList_New(dataset->n_atoms);
  for (i = 0; i < dataset->n_atoms; i++) {
    PyList_SetItem(mapping_to_primitive, i,
                   PyLong_FromLong((long) dataset->mapping_to_primitive[i]));
  }
  
  PyList_SetItem(array, n, mapping_to_primitive);
  n++;


  /* Standardized unit cell */
  primitive_types = PyList_New(dataset->n_std_atoms);
  positions = PyList_New(dataset->n_std_atoms);
  magmom = PyList_New(dataset->n_std_atoms);
  for (i = 0; i < dataset->n_std_atoms; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->positions[i][j]));
    }
    PyList_SetItem(positions, i, vec);

    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->magmom[i][j]));
    }
    PyList_SetItem(magmom, i, vec);

    PyList_SetItem(primitive_types, i, PyLong_FromLong((long) dataset->primitive_types[i]));
  }
  PyList_SetItem(array, n, primitive_types);
  n++;
  PyList_SetItem(array, n, positions);
  n++;
  PyList_SetItem(array, n, magmom);
  n++;


  /* Point group */
  /* PyList_SetItem(array, n, PyLong_FromLong((long) dataset->pointgroup_number)); */
  /* n++; */
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->pointgroup_symbol));
  n++;

  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->msgtype));
  n++;

  assert(n == len_list);

  spg_free_magneticdataset(dataset);

  return array;
}



static PyObject * py_get_error_message(PyObject *self, PyObject *args)
{
  MspglibError error;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  error = mspg_get_error_code();

  return PYUNICODE_FROMSTRING(mspg_get_error_message(error));
}
