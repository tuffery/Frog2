# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _vflib
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class ARGEdit(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ARGEdit, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ARGEdit, name)
    def InsertNode(*args): return apply(_vflib.ARGEdit_InsertNode,args)
    def DeleteNode(*args): return apply(_vflib.ARGEdit_DeleteNode,args)
    def DeleteEdge(*args): return apply(_vflib.ARGEdit_DeleteEdge,args)
    def InsertEdge(*args): return apply(_vflib.ARGEdit_InsertEdge,args)
    def NodeCount(*args): return apply(_vflib.ARGEdit_NodeCount,args)
    def __init__(self,*args):
        self.this = apply(_vflib.new_ARGEdit,args)
        self.thisown = 1
    def __del__(self, destroy= _vflib.delete_ARGEdit):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C ARGEdit instance at %s>" % (self.this,)

class ARGEditPtr(ARGEdit):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = ARGEdit
_vflib.ARGEdit_swigregister(ARGEditPtr)

class GraphMatcher(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GraphMatcher, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GraphMatcher, name)
    def __init__(self,*args):
        self.this = apply(_vflib.new_GraphMatcher,args)
        self.thisown = 1
    def matchVF(*args): return apply(_vflib.GraphMatcher_matchVF,args)
    def matchVFMono(*args): return apply(_vflib.GraphMatcher_matchVFMono,args)
    def matchVFSub(*args): return apply(_vflib.GraphMatcher_matchVFSub,args)
    def matchVF2(*args): return apply(_vflib.GraphMatcher_matchVF2,args)
    def matchVF2Sub(*args): return apply(_vflib.GraphMatcher_matchVF2Sub,args)
    def matchVF2Mono(*args): return apply(_vflib.GraphMatcher_matchVF2Mono,args)
    def matchUll(*args): return apply(_vflib.GraphMatcher_matchUll,args)
    def matchUllSub(*args): return apply(_vflib.GraphMatcher_matchUllSub,args)
    def __del__(self, destroy= _vflib.delete_GraphMatcher):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C GraphMatcher instance at %s>" % (self.this,)

class GraphMatcherPtr(GraphMatcher):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = GraphMatcher
_vflib.GraphMatcher_swigregister(GraphMatcherPtr)


