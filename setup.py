# coding=utf-8

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

from shutil import which

import io
import os
import re
import sys
import shlex
import subprocess

# CMake arguments
cmake_defines = "-DVIPSTER_PYLIB=YES"
if os.name == 'nt':
    # force msys makefile generation
    cmake_defines += " -G \"MSYS Makefiles\""
else:
    # force static linking
    cmake_defines += " -DBUILD_SHARED_LIBS=OFF"
# circumvent hypot definition error in 3.x < 3.7
if sys.version_info.major == 3 and sys.version_info.minor < 7:
    cmake_defines += " -DCMAKE_CXX_FLAGS=\"-D_hypot=hypot\""


# find CMake
CMAKE_EXE = os.environ.get('CMAKE_EXE', which('cmake'))
if not CMAKE_EXE:
    print('CMake executable not found.'
          'Set CMAKE_EXE environment or update your PATH.')
    sys.exit(1)

# on windows, find gcc-libraries as long as we cannot build with msvc
if os.name == 'nt':
    extralibs = [
        ('../../vipster',
         [t for t in map(which,
          ['libwinpthread-1.dll',
           'libgcc_s_seh-1.dll',
           'libgcc_s_dw2-1.dll',
           'libstdc++-6.dll',
           ]) if t])]
else:
    extralibs = []


# CMake* classes borrowed from https://github.com/raydouglass/cmake_setuptools
# until this works with other generators/python2
class CMakeExtension(Extension):
    """setuptools.Extension for cmake."""

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuildExt(build_ext):
    """setuptools build_ext which builds using cmake & make."""

    def build_extension(self, ext):
        if isinstance(ext, CMakeExtension):
            output_dir = os.path.abspath(
                os.path.dirname(self.get_ext_fullpath(ext.name)))

            build_type = 'Debug' if self.debug else 'Release'
            cmake_args = [CMAKE_EXE,
                          ext.sourcedir,
                          '-DBUILD_TESTING=OFF',
                          '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + output_dir,
                          '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + output_dir,
                          '-DCMAKE_BUILD_TYPE=' + build_type]
            cmake_args.extend(
                [x for x in
                 shlex.split(cmake_defines)
                 if x])

            env = os.environ.copy()
            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)
            subprocess.check_call(cmake_args,
                                  cwd=self.build_temp,
                                  env=env)
            subprocess.check_call([CMAKE_EXE, '--build', '.', '-j',
                                   '--target', ext.name],
                                  cwd=self.build_temp,
                                  env=env)
            print()
        else:
            super().build_extension(ext)


# parse and prepare readme
def readfile(arg):
    here = os.path.abspath(os.path.dirname(__file__))
    with io.open(os.path.join(here, arg)) as f:
        return f.read()


readme = readfile('README.md')
readme = re.sub('(INSTALL.md)', r'https://github.com/sgsaenger/'
                r'vipster/blob/master/\1', readme)
readme = re.sub('(util/vipster.png)', r'https://raw.githubusercontent.com/'
                r'sgsaenger/vipster/master/\1', readme)

setup(name="vipster",
      version=re.findall(r'project\(Vipster VERSION ([0-9.]*)',
                         readfile('CMakeLists.txt'))[0],
      author="Sebastian Gsänger",
      url="https://github.com/sgsaenger/vipster",
      description="A pre- and post-processing toolkit "
                  "for atomistic simulations.",
      long_description=readme,
      long_description_content_type="text/markdown",
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Programming Language :: C++",
          "Programming Language :: Python :: Implementation :: CPython",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Scientific/Engineering :: Physics",
      ],
      python_requires='>=3.7',
      keywords=['chemistry'],
      license="GPL",
      ext_modules=[CMakeExtension('pyvipster', '.')],
      data_files=extralibs,
      cmdclass={'build_ext': CMakeBuildExt},
      zip_safe=False,
      )
