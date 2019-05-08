# coding=utf-8

from setuptools import setup
from cmake_setuptools import CMakeExtension, CMakeBuildExt

import io
from os import path, environ, name as osname
import re

environ["CMAKE_COMMON_VARIABLES"] = "-DPYTHON=YES -DPYPI=YES"
if osname == 'nt':
    # force mingw makefile generation
    environ["CMAKE_COMMON_VARIABLES"] += " -G \"MinGW Makefiles\""
else:
    # force static linking
    environ["CMAKE_COMMON_VARIABLES"] += " -DBUILD_SHARED_LIBS=NO"
here = path.abspath(path.dirname(__file__))


def readfile(arg):
    with io.open(path.join(here, arg)) as f:
        return f.read()


readme = readfile('README.md')
readme = re.sub('(INSTALL.md)', r'https://github.com/sgsaenger/vipster/blob/master/\1', readme)
readme = re.sub('(dist/vipster.png)', r'https://raw.githubusercontent.com/sgsaenger/vipster/master/\1', readme)

setup(
        name="vipster",
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
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            ],
        keywords=['chemistry'],
        license="GPL",
        ext_modules=[CMakeExtension('all', '.')],
        cmdclass={'build_ext': CMakeBuildExt},
        )
