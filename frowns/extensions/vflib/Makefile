TOP        = ../..
SWIG       = $(TOP)/../swig
CXXSRCS    = vflib.cxx vflib-2.0/src/argedit.cxx vflib-2.0/src/argloader.cxx vflib-2.0/src/argraph.cxx vflib-2.0/src/error.cxx vflib-2.0/src/gene.cxx vflib-2.0/src/gene_mesh.cxx vflib-2.0/src/match.cxx vflib-2.0/src/sd_state.cxx vflib-2.0/src/sortnodes.cxx vflib-2.0/src/ull_state.cxx vflib-2.0/src/ull_sub_state.cxx vflib-2.0/src/vf2_mono_state.cxx vflib-2.0/src/vf2_state.cxx vflib-2.0/src/vf2_sub_state.cxx vflib-2.0/src/vf_mono_state.cxx vflib-2.0/src/vf_state.cxx vflib-2.0/src/vf_sub_state.cxx vflib-2.0/src/xsubgraph.cxx
OBJS       = vflib.o argedit.o argloader.o argraph.o error.o gene.o gene_mesh.o match.o sd_state.o sortnodes.o ull_state.o ull_sub_state.o vf2_mono_state.o vf2_state.o vf2_sub_state.o vf_mono_state.o vf_state.o vf_sub_state.o xsubgraph.o
INCLUDES   = vflib-2.0/include
TARGET     = vflib
INTERFACE  = vflib.i
LIBS       = -lm

all::
	$(MAKE) -f $(TOP)/Makefile CXXSRCS='$(CXXSRCS)' SWIG='$(SWIG)' \
	INCLUDE='-Ivflib-2.0/include' OBJS='$(OBJS)' TARGET='$(TARGET)' INTERFACE='$(INTERFACE)' python_cpp

static::
	$(MAKE) -f $(TOP)/Makefile CXXSRCS='$(CXXSRCS)' SWIG='$(SWIG)' \
	TARGET='mypython' INTERFACE='$(INTERFACE)' python_cpp_static

clean::
	$(MAKE) -f $(TOP)/Makefile python_clean
	rm -f vflib.py

check: all
