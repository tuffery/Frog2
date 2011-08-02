#include <iostream>
#include <string>
#include <list>

#include "Python.h"
#include "pysssr.h"

#define INPUTERROR "Input not the number of atoms followed a list of pairs of integers [(1,2), (3,4)...]"

#define THROWINPUTERROR PyErr_SetString(PyExc_TypeError, INPUTERROR); return NULL;
#define THROW_RUNTIME(x) PyErr_SetString(PyExc_RuntimeError, x);  return NULL;

Molecule * get_mol() {
  // OC[CH]1O[CH](SC2=C(Cl)C=CC(=C2)Cl)[CH](O)[CH](O)[CH]1O
  const int N=20;
  Molecule *mol = new Molecule(N);
  
  
  mol->add_bond(0, 1);
  mol->add_bond(1, 2);
  mol->add_bond(2, 3);
  mol->add_bond(3, 4);
  mol->add_bond(4, 5);
  mol->add_bond(5, 6);
  mol->add_bond(6, 7);
  mol->add_bond(7, 8);
  mol->add_bond(7, 9);
  mol->add_bond(9, 10);
  mol->add_bond(10, 11);
  mol->add_bond(11, 12);
  mol->add_bond(6, 12);
  mol->add_bond(11, 13);
  mol->add_bond(4, 14);
  mol->add_bond(14, 15);
  mol->add_bond(14, 16);
  mol->add_bond(16, 17);
  mol->add_bond(16, 18);
  mol->add_bond(2, 18);
  mol->add_bond(18, 19);
  
  return mol;
}

PyObject *get_rings(RingSet &ringSet) {
  //assert(0);
  // For speed purposes, this also deletes the ringSet as it
  // consumes it.
  PyObject *rings = PyTuple_New(ringSet.size());
  std::list<Ring >::iterator iter;
  int ringIndex = 0;
  
  for(iter = ringSet.rings.begin(); iter != ringSet.rings.end(); iter++, ringIndex++) {
    Ring &ring = *iter;
    PyObject *verts = PyTuple_New(ring.count());
    PyTuple_SetItem(rings, ringIndex, verts);
    
    unsigned int lastvert = 0;
    //    std::cout << "number of atoms " << ringSet.mol->numberOfAtoms << "\n";
    for(int index=0; index < ringSet.mol->numberOfAtoms; index++) {
      if(ring[index]) 
	PyTuple_SetItem(verts, lastvert++, PyInt_FromLong(index));
    }

    if (lastvert != ring.count()) {
      std::cout << "\t" << ring << std::endl;
      std::cout << "vertices and expected vertices don't match\n";
      std::cout << "lastvert " << lastvert << " ring count " << ring.count() << "\n";

      THROW_RUNTIME("Number of computed vertices and number of expected vertices doesn't match.");
    }

  }
  return rings;
}

PyObject *sssr(PyObject *numAtoms, PyObject *atoms) {
  if (!PyInt_Check(numAtoms)) {
    THROWINPUTERROR;
  }
  
  int natoms = PyInt_AS_LONG(numAtoms);
  
  Molecule *mol = new Molecule(natoms);
  
  if (PySequence_Check(atoms) == 0) {
    // Oops not a sequence, so throw an error
    THROWINPUTERROR;
  }
  
  // go throw all the sequences
  PyObject *bond, *from, *to;
  int atom1, atom2;
  int atomLength = PySequence_Length(atoms);

  // If we have more atoms than MAXATOMS then throw
  // a runtime error
  if (atomLength > MAXATOMS) {
    THROW_RUNTIME("pysssr can only handle this many atoms.");
  }

  for (int i=0; i<PySequence_Length(atoms); i++) {
    bond = PySequence_GetItem(atoms, i);
    
    if (PySequence_Check(bond) == 0 || PySequence_Length(bond) != 2) {
      THROWINPUTERROR;
    }
    
    from = PySequence_GetItem(bond, 0);
    to = PySequence_GetItem(bond, 1);
    
    // If either of the from or to atoms is not an integer, this is an input
    // error
    if (!PyInt_Check(from) || !PyInt_Check(to)) {
      THROWINPUTERROR;
    }
    
    atom1 = PyInt_AS_LONG(from);
    atom2 = PyInt_AS_LONG(to);
    
    Py_DECREF(bond);
    Py_DECREF(from);
    Py_DECREF(to);
    
    if(atom1 > natoms || atom1 < 0 || atom2 > natoms || atom2 < 0) {
      //      std::cout << "Type Error should be here\n";
      PyErr_SetString(PyExc_TypeError, 
		      "atom index out of range of number of atoms"); 
      return Py_None;
    }
    
    mol->add_bond(atom1, atom2);
    
  }
  
  
  try {
    RingSet &result = ringDetection(mol, -1);  
    PyObject *rings = get_rings(result);
    
    delete &result;
    delete mol;
    return rings;
  } catch (...) {
    THROW_RUNTIME("Ring detection failed\n");
  }


}

static PyObject *_pysssr(PyObject *self, PyObject *args)
{
  PyObject *numAtoms, *atoms;

  //  std::cout << "parsing tuples\n";
  if (!PyArg_ParseTuple(args, "OO", &numAtoms, &atoms))
    return NULL;

  //  std::cout << "got " << PyInt_AS_LONG(numAtoms) << "atoms\n";
  return sssr(numAtoms, atoms);
}


static PyMethodDef SSSRMethods[] = {
    {"sssr",  _pysssr, METH_VARARGS, "given a list of atoms and connections, return the smallest set of smallest rings."},
    {NULL, NULL}        /* Sentinel */
};

extern "C" {
  void
  init_pysssr(void)
  {
    //    std::cout << "Loading _pysssr version 1.0\n";
    (void) Py_InitModule("_pysssr", SSSRMethods);
  }

}
