#ifndef __VFLIB__
#define __VFLIB__

#include <iostream>
#include <string>
#include <list>

#include "error.h"
#include "Python.h"

// Basic graph structure files
#include "argraph.h"
#include "argedit.h"
#include "match.h"

// Isomorphism files
#include "vf_state.h"
#include "vf_sub_state.h"
#include "vf_mono_state.h"

#include "vf2_state.h"
#include "vf2_sub_state.h"
#include "vf2_mono_state.h"

#include "ull_state.h"
#include "ull_sub_state.h"

const int DEBUG=0;  
typedef ARGraph<PyObject, PyObject> PyGraph;

/////////////////////////////////////////////////////////////////////
//  PyARGEdit - python wrapper class for ARGEdit allows
//               node attributes and edge attributes to be
//               python objects
/////////////////////////////////////////////////////////////////////
struct edge {
  // Just a little data structure to store edge information.
  int src;
  int dst;
  edge(int a, int b) { src=a; dst=b; }
};

class PyARGEdit:public ARGEdit {
public:
  int pyargedit_num_edges;
  std::list<edge> edges;               // Edges are kept because the isomorphism
  //  programs return the mapped nodes.
  //  we will want to return the mapped
  //  edges as well...
  // WARNING!  There can only be two edge
  //  between a two nodes edge(a,b) and edge(b,a)
  //  multiple links between nodes are not supported!
  std::list<PyObject *> pyvertices;    // We need to keep references to the
  std::list<PyObject *> pyedges;       //  edges and vertices so they don't
  //  go away under our butts.
  
  PyARGEdit():ARGEdit() {
    pyargedit_num_edges = 0;
  }

  ~PyARGEdit() {
    // Clean up the saved vertices and edges
    if (DEBUG) printf("deleting pyargedit\n");
    std::list<PyObject *>::iterator i;// = pyvertices.begin();
    
    // We can decref the saved vertices and edges here to
    // clean up memory
    for(i=pyvertices.begin(); i!=pyvertices.end(); i++) {
      Py_DECREF(*i);
      *i=NULL;
    }
    
    for(i=pyedges.begin(); i!=pyedges.end(); i++) {
      Py_DECREF(*i);
      *i=NULL;
    }
  }
  
  node_id PyInsertNode(PyObject *ob);
  void PyDeleteNode(node_id id);
  void PyDeleteEdge(node_id id1, node_id id2);
  void PyInsertEdge(node_id id1, node_id id2, PyObject *ob);
  int PyNodeCount();
};

class PyObjectCompareError {};
bool pyobject_compare(PyObject *a, PyObject *b);
bool collectionWrapper(int n, node_id ni1[], node_id ni2[], void *usr_data);

///////////////////////////////////////////////////////////
// A collection is a data object that holds the results
// of a match in a python list.  The only reason that
// this object exists is to dereference the reult
// object due to some wierdness in boost.  When a
// PyObject is returned from a boost call, the reference
// automatically increased by one.
// However, using PyList_New starts with a reference count of
// one!  This means that returning anything created from PyList_New
// via boost will never be garbage collected.
//  So we need a mechanism to clean this up.  When this class
//  is deleted by C++, the list reference count is reduced by one.
//  So this is a pretty darn big hack.
//
// You can also set the limit of results here.
// atLimit() returns 1 if the limit has been reached.
// next() increments the current count.
class GraphMatcher;

class Collection {
public:
  int limit;
  int count;
  PyObject *result;
  GraphMatcher &matcher;
  
  Collection(int maxCount, GraphMatcher &callingGraphMatcher):matcher(callingGraphMatcher) {
    limit = maxCount;
    count = 0;
    result = PyList_New(0);
  }
  
  Collection::~Collection() {
    //Py_DECREF(result);
  }
  
  inline bool atLimit() {
    if(count == -1) return false;
    return count == limit;
  }
  
  inline void next() { count++; }
};

////////////////////////////////////////////////////////////
// a wrapper for the graph matchers
// usage will be
// matcher = GraphMatcher(g) # g is a PyARGEdit object
// matcher.matchVF2(h) # h is a PyARGEdit object -> 
//   returns a list of ((node), (edge)) mappings for each match found


class GraphMatcher {
public:
  PyGraph *g, *_h;
  PyARGEdit *original_match;
  Collection *_collection;
  // The following need to be copied from incoming match
  //  object because python may try to reclaim the instance
  std::list<edge> edges;																																																																				
  std::list<PyObject *> pyvertices;		
  std::list<PyObject *> pyedges;		
  
  GraphMatcher(PyARGEdit &match) {
    std::list<PyObject *>::iterator i;
    
    g = new PyGraph(&match);
    g->SetNodeCompat(pyobject_compare);
    g->SetEdgeCompat(pyobject_compare);
    
    original_match = &match;
    // INREF the nodes and edges
    for(i=match.pyvertices.begin(); i!= match.pyvertices.end(); i++) {
      pyvertices.push_back(*i);
      Py_INCREF(*i);
    }
    
    for(i=match.pyedges.begin(); i!= match.pyedges.end(); i++) {
      pyvertices.push_back(*i);
      Py_INCREF(*i);
    }
    
    std::list<edge>::iterator edge_iterator;
    for(edge_iterator=match.edges.begin(); edge_iterator!=match.edges.end(); edge_iterator++) {
      edges.push_back(edge((*edge_iterator).src, (*edge_iterator).dst));
    }
    
    _h = NULL;
    _collection = NULL;
  }
  
  ~GraphMatcher() {
    std::list<PyObject *>::iterator i;
    
    if (DEBUG) printf("deleting graph matcher\n");
    delete g;
    if(_h) delete _h;
    if(_collection) delete _collection;
    
    // DECREF the nodes and edges for the match object
    for(i=pyvertices.begin(); i!= pyvertices.end(); i++)
      Py_DECREF(*i);
    
    for(i=pyedges.begin(); i!= pyedges.end(); i++)
      Py_DECREF(*i);
    
  }
  
  inline void prepareMatch(PyARGEdit &target) {	  
    if(_h) delete _h;
    _h = new PyGraph(&target);
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
  PyObject *matchVF(PyARGEdit &target, int limit);
  PyObject *matchVFMono(PyARGEdit &target, int limit);
  PyObject *matchVFSub(PyARGEdit &target, int limit);
  
  // VF2 mathcers (isomorphic, subgraph isomorphic and monoisomorphic)
  PyObject *matchVF2(PyARGEdit &target, int limit);
  PyObject *matchVF2Sub(PyARGEdit &target, int limit);
  PyObject *matchVF2Mono(PyARGEdit &target, int limit);
  
  // Ullman matchers (Isomorphic and subgraph isomorphic)
  PyObject *matchUll(PyARGEdit &target, int limit);  
  PyObject *matchUllSub(PyARGEdit &target, int limit);
};
#endif
