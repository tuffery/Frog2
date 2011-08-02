/* file: vflib.i */
%module vflib

%{
#include "vflib.h"
  %}

typedef int node_id;

%name(ARGEdit) class PyARGEdit {
public:
  %name(InsertNode) node_id PyInsertNode(PyObject *ob);
  %name(DeleteNode) void PyDeleteNode(node_id id);
  %name(DeleteEdge) void PyDeleteEdge(node_id id1, node_id id2);
  %name(InsertEdge) void PyInsertEdge(node_id id1, node_id id2, PyObject *ob);
  %name(NodeCount) int PyNodeCount();
};

class GraphMatcher {
 public:
  GraphMatcher(PyARGEdit &match);
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
