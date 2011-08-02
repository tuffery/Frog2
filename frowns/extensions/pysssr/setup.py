import distutils, os, sys
from distutils.core import setup, Extension
from distutils.command import build_ext

src = [ 'pysssr.cxx',
        'sssr.cxx'
        ]

objs = [
  'pysssr.cxx',
  'sssr.cxx'
  ]


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

setup(name = "_pysssr",
      version = "1.0",
      cmdclass = {'build_ext': my_build_ext, },
      ext_modules = [Extension("_pysssr", src)])




