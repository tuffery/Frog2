"""Python setup file for vflib

SWIG 1.3+ must be installed for this to build correctly.
"""
import distutils, os, sys
from distutils.core import setup, Extension
from distutils.command import build_ext

src = ['vflib.cxx', 'vflib-2.0/src/argedit.cxx', 'vflib-2.0/src/argloader.cxx',
       'vflib-2.0/src/argraph.cxx', 'vflib-2.0/src/error.cxx', 'vflib-2.0/src/gene.cxx',
       'vflib-2.0/src/gene_mesh.cxx', 'vflib-2.0/src/match.cxx',
       'vflib-2.0/src/sd_state.cxx', 'vflib-2.0/src/sortnodes.cxx',
       'vflib-2.0/src/ull_state.cxx', 'vflib-2.0/src/ull_sub_state.cxx',
       'vflib-2.0/src/vf2_mono_state.cxx', 'vflib-2.0/src/vf2_state.cxx',
       'vflib-2.0/src/vf2_sub_state.cxx', 'vflib-2.0/src/vf_mono_state.cxx',
       'vflib-2.0/src/vf_state.cxx', 'vflib-2.0/src/vf_sub_state.cxx',
       'vflib-2.0/src/xsubgraph.cxx']

objs = ['vflib.o', 'argedit.o', 'argloader.o', 'argraph.o', 'error.o', 'gene.o',
        'gene_mesh.o', 'match.o', 'sd_state.o', 'sortnodes.o', 'ull_state.o',
        'ull_sub_state.o', 'vf2_mono_state.o', 'vf2_state.o', 'vf2_sub_state.o',
        'vf_mono_state.o', 'vf_state.o', 'vf_sub_state.o', 'xsubgraph.o']


class my_build_ext(build_ext.build_ext):
   def swig_sources (self, sources):

        """Walk the list of source files in 'sources', looking for SWIG
        interface (.i) files.  Run SWIG on all that are found, and
        return a modified 'sources' list with SWIG source files replaced
        by the generated C (or C++) files.
        """
        self.swig_cpp = 1
        new_sources = []
        swig_sources = []
        swig_targets = {}

        # XXX this drops generated C/C++ files into the source tree, which
        # is fine for developers who want to distribute the generated
        # source -- but there should be an option to put SWIG output in
        # the temp dir.

        if self.swig_cpp:
            target_ext = '.cxx'
        else:
            target_ext = '.c'

        for source in sources:
            (base, ext) = os.path.splitext(source)
            if ext == ".i":             # SWIG interface file
                # new_sources.append(base + target_ext)
                new_sources.append(base + '_wrap' + target_ext)  #new
                
                swig_sources.append(source)
                swig_targets[source] = new_sources[-1]
            else:
                new_sources.append(source)

        if not swig_sources:
            return new_sources

        swig = self.find_swig()
        # swig_cmd = [swig, "-python", "-dnone", "-ISWIG"]
        swig_cmd = [swig, "-python"]  #new
        
        if self.swig_cpp:
            swig_cmd.append("-c++")

        for source in swig_sources:
            target = swig_targets[source]
            self.announce("swigging %s to %s" % (source, target))
            self.spawn(swig_cmd + ["-o", target, source])

        return new_sources

    # swig_sources ()  

if sys.platform == "win32":
  ########   WINDOWS  ###########################################
  extra_compile_args = ["/GR"]
else:
  ########  UNIX      ###########################################
  ## fix for boost and g++
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

setup(name = "PySomorph",
      version = "1.1",
      py_modules=['vflib'],
      cmdclass = {'build_ext': my_build_ext, },
      ext_modules = [Extension("_vflib", ["vflib.i"]+src)], include_dirs=["vflib-2.0/include"])




