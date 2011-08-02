///////////////////////////////////////////////////////////////////////
// vflib version 1.0
// See the enclosed license for use (vflib.txt)
// Wrapper created by Brian Kelley using boost python 
//  http://www.boost.org
///////////////////////////////////////////////////////////////////////
#include "vflib.h"

void memtest() {
  // This just loops over graph creation to ensure that memory
  // is not being leaked.
  ARGEdit *g;
  Graph *graph;
  while(1) {
    // create a ArgEdit structure
    g = new ARGEdit();
    graph = new Graph(g);
    printf("%x", graph->NodeCount());
    delete g;
    delete graph;
  }
}

node_id PyARGEdit::PyInsertNode(PyObject *ob) { 
  Py_INCREF(ob); // incref the object so it won't go away beneath us.
  pyvertices.push_back(ob);
  return InsertNode((void *) ob); 
}

void PyARGEdit::PyDeleteNode(node_id id) {
  // Danger, this might throw a run time error
  if (id > pyvertices.size()) {
    // raise 
  }
  DeleteNode(id);
}

void PyARGEdit::PyDeleteEdge(node_id id1, node_id id2) {
  // Danger, this might throw a run time error
  DeleteEdge(id1, id2);
}

void PyARGEdit::PyInsertEdge(node_id id1, node_id id2, PyObject *ob) {
  Py_INCREF(ob); // incref the object so it won't go away beneath us.
  pyedges.push_back(ob);
  pyargedit_num_edges++;
  edges.push_back(edge(id1,id2));
  InsertEdge(id1, id2, (void *) ob);
}

int PyARGEdit::PyNodeCount() {
    return NodeCount();
}

////////////////////////////////////////////////////////////
// Test for returning a newly created by object
PyObject *getObject() {
  PyObject *t = PyTuple_New(1);
  PyTuple_SetItem(t, 0, PyInt_FromLong(1));
  PyObject *result = PyList_New(1);
  PyList_SetItem(result, 0, t);
  return result;
}
////////////////////////////////////////////////////////////
// a comparison function for python based node and edge
// attributes
bool pyobject_compare(PyObject *a, PyObject *b) {
  // return 1 if the objects are compareable
  // return 0 if not and let Python do all the work!!!
  
  PyObject *error;
  
  int result = PyObject_Compare(a, b);
  error = PyErr_Occurred();
  if(error) {
    // if there is a python error while executing PyObject_Compare
    //  we can trap it and send it back to the interpreter
    // This code actually sends a PyCompare error becuase I'm not
    // quite sure how to raise the original error yet.
    //  however, I do print out the original error so the user can
    //  tell where it's happened.
    PyErr_Print();
    PyErr_SetObject(PyExc_RuntimeError, PyString_FromString("PyCompare error"));//BOOST_PYTHON_CONVERSION::to_python("PyCompare error"));				 return false;
    throw PyObjectCompareError();
  }
  
  if (!result) {
    return true;
  } else {
    return false;
  }	
}

// This is the match visitor defined by the VFLIB sources
// n is the number of edges, n1 and n2 are vertices in query->match
bool collectionWrapper(int n, node_id ni1[], node_id ni2[], void *usr_data) {
  // PyList *result = (PyList *)usr_data;
  Collection *result = (Collection *)usr_data;
  std::list<edge>::iterator i, end;
  unsigned int index;
  // collect the nodes of the result

  // Don't collect empty results
  if (!n) return false;

  PyObject *nodes = PyTuple_New(n);
  
  for(index=0; index<(unsigned int)n; index++) {
    int gposition = ni1[index];
    int hposition = ni2[index];
    assert(gposition<n);
    PyObject *o = result->matcher._h->GetNodeAttr(hposition);
    Py_INCREF(o);	  
    PyTuple_SetItem(nodes, gposition, o);
  }
  
  // collect the edges of the result
  
  PyObject *edges = PyTuple_New(result->matcher.edges.size());
  
  end = result->matcher.edges.end();
  for(i = result->matcher.edges.begin(), index = 0; i != end; i++, index++) {			  
    edge &e = *i;			  
    node_id from = e.src;			  
    node_id to = e.dst;			  
    
    if (!result->matcher._h->HasEdge(ni2[from], ni2[to], NULL)){
      //PyErr_SetObject(PyExc_AssertionError, BOOST_PYTHON_CONVERSION::to_python(index));				  
      //throw boost::python::error_already_set();
    }
    
    PyObject *attribute = result->matcher._h->GetEdgeAttr(ni2[from], ni2[to]);
    
    assert(index<result->matcher.edges.size());
    Py_INCREF(attribute);
    PyTuple_SetItem(edges, index, attribute);
  }
  
  
  PyObject *nodes_and_edges = PyTuple_New(2);
  
  PyTuple_SetItem(nodes_and_edges, 0, nodes);
  PyTuple_SetItem(nodes_and_edges, 1, edges);
  
  
  PyList_Append(result->result, nodes_and_edges);
  Py_DECREF(nodes_and_edges);
  
  
  result->next();
  
  if(result->atLimit()) return true;
  return false;
}


// Matchers work in the following way:
// first prepare the target to be matched
// create the state object that actually does the matching
// create a collection object to show the results (this is
//  the data object sent when a match is found, see 
//  the collectionWrapper function below)
// return the _collection->result
//
// There's a lot a bit of replicated code but what the hey!

// VF2 matchers (isomorphic, subgraph isomorphic and monoisomorphic)
PyObject *GraphMatcher::matchVF(PyARGEdit &target, int limit) {
  prepareMatch(target);
  VFState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}
  
PyObject *GraphMatcher::matchVFMono(PyARGEdit &target, int limit) {
  prepareMatch(target);
  VFMonoState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}

PyObject *GraphMatcher::matchVFSub(PyARGEdit &target, int limit) {
  prepareMatch(target);
  VFSubState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}

// VF2 mathcers (isomorphic, subgraph isomorphic and monoisomorphic)
PyObject *GraphMatcher::matchVF2(PyARGEdit &target, int limit) {
  prepareMatch(target);
  VF2State s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}

PyObject *GraphMatcher::matchVF2Sub(PyARGEdit &target, int limit) {
  prepareMatch(target);
  VF2SubState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}

PyObject *GraphMatcher::matchVF2Mono(PyARGEdit &target, int limit) {
  prepareMatch(target);
  VF2MonoState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}


// Ullman matchers (Isomorphic and subgraph isomorphic)
PyObject *GraphMatcher::matchUll(PyARGEdit &target, int limit) {
  prepareMatch(target);
  UllState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}

PyObject *GraphMatcher::matchUllSub(PyARGEdit &target, int limit) {
  prepareMatch(target);
  UllSubState s0(g, _h);
  if(_collection) delete _collection;
  
  _collection = new Collection(limit, *this);
  
  match(&s0, collectionWrapper, _collection);
  return _collection->result;
}	


// BOOST_PYTHON_MODULE_INIT(vflib)
// {
//   try
//   {
//     /*******************************************************/
//     // Test code
//     // Create an object representing this extension module.
//     python::module_builder this_module("vflib");

//     this_module.def(memtest, "memtest");
// 	this_module.def(getObject, "getObject");

// 	python::class_builder<PyARGEdit> argedit_class(this_module, "ARGEdit");
// 	argedit_class.def(python::constructor<>());
// 	argedit_class.def(&PyARGEdit::PyInsertNode, "InsertNode");
// 	argedit_class.def(&PyARGEdit::PyInsertEdge, "InsertEdge");
// 	argedit_class.def(&PyARGEdit::PyDeleteNode, "DeleteNode");
// 	argedit_class.def(&PyARGEdit::PyDeleteEdge, "DeleteEdge");
// 	argedit_class.def(&PyARGEdit::PyNodeCount, "NodeCount");


// 	python::class_builder<GraphMatcher> graphmatcher_class(this_module, "GraphMatcher");
// 	graphmatcher_class.def(python::constructor<PyARGEdit &>());

// 	graphmatcher_class.def(&GraphMatcher::matchVF, "matchVF");
// 	graphmatcher_class.def(&GraphMatcher::matchVFMono, "matchVFMono");
// 	graphmatcher_class.def(&GraphMatcher::matchVFSub, "matchVFSub");

// 	graphmatcher_class.def(&GraphMatcher::matchVF2, "matchVF2");
// 	graphmatcher_class.def(&GraphMatcher::matchVF2Mono, "matchVF2Mono");
// 	graphmatcher_class.def(&GraphMatcher::matchVF2Sub, "matchVF2Sub");

// 	graphmatcher_class.def(&GraphMatcher::matchUll, "matchUll");
// 	graphmatcher_class.def(&GraphMatcher::matchUllSub, "matchUllSub");
	
// 	  }
//   catch(...)
//   {
//     python::handle_exception(); // Deal with the exception for Python
//   }
// }
