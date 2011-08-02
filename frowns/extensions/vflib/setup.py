"""Python setup file for vflib

See setup_swig.py if you want to install using swig.
"""
import distutils, os, sys
from distutils.core import setup, Extension


if sys.platform == "win32":
  ########   WINDOWS  ###########################################
  pass
else:
  ########  UNIX      ###########################################
  ## fix for boost and g++
  extra_compile_args = ["-O3"]
  from distutils import sysconfig
  save_init_posix = sysconfig._init_posix
  def my_init_posix():
    print "my_init_posix: changing gcc to g++"
    save_init_posix()
    g = sysconfig._config_vars
    g['CC'] = 'g++'
    if sys.platform == "darwin":
      g['LDSHARED'] = 'g++ -bundle -undefined suppress -flat_namespace'
    else:
      g['LDSHARED'] = 'g++ -shared'
  sysconfig._init_posix = my_init_posix

src = ['vflib.cxx', 'vflib_wrap.cxx',
       'vflib-2.0/src/argedit.cxx', 'vflib-2.0/src/argloader.cxx',
       'vflib-2.0/src/argraph.cxx', 'vflib-2.0/src/error.cxx', 'vflib-2.0/src/gene.cxx',
       'vflib-2.0/src/gene_mesh.cxx', 'vflib-2.0/src/match.cxx',
       'vflib-2.0/src/sd_state.cxx', 'vflib-2.0/src/sortnodes.cxx',
       'vflib-2.0/src/ull_state.cxx', 'vflib-2.0/src/ull_sub_state.cxx',
       'vflib-2.0/src/vf2_mono_state.cxx', 'vflib-2.0/src/vf2_state.cxx',
       'vflib-2.0/src/vf2_sub_state.cxx', 'vflib-2.0/src/vf_mono_state.cxx',
       'vflib-2.0/src/vf_state.cxx', 'vflib-2.0/src/vf_sub_state.cxx',
       'vflib-2.0/src/xsubgraph.cxx']

setup(name = "PySomorph",
      version = "1.1",
      py_modules=['vflib'],
      ext_modules = [Extension("_vflib", src)],
      include_dirs=["vflib-2.0/include"])




