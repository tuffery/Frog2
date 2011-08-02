from distutils.core import setup, Extension
import sys, os

setup(name="frowns",
      version="0.9a",
      author="Brian Kelley",
      author_email = "bkelley@wi.mit.edu",
      url="http://frowns.sourceforge.net/",
      package_dir = {'':".."},
      packages=['frowns',
                'frowns/Canonicalization',
                'frowns/Chirality',
                'frowns/Depict',
                'frowns/Fingerprint',
                'frowns/Objects',
                'frowns/perception',
                'frowns/mdl_parsers',
                'frowns/smiles_parsers',
                'frowns/Smarts',
                'frowns/Utils',
                ]
      )



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

# This builds the pysssr (smallest set of smallest rings) extension
os.chdir("extensions/pysssr")
src = [ 'pysssr.cxx', 'sssr.cxx' ]
    
setup(name = "_pysssr",
      version = "1.1",
      ext_modules = [Extension("_pysssr", src)])

os.chdir("../vflib")

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


