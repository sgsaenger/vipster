#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy as np

mol_c = Extension(name="vipster.mol_c",
                  sources=["vipster/mol_c.c"],
                  include_dirs=[np.get_include()],
                  extra_compile_args=['-std=c99'])
gui_c = Extension(name="vipster.gui.gui_c",
                  sources=["vipster/gui/gui_c.c"],
                  include_dirs=[np.get_include()],
                  extra_compile_args=['-std=c99'])

setup(
    name="vipster",
    version="0.8.0",
    description="VIsual Periodic STructure EditoR",
    author="Sebastian Gsänger",
    author_email="sebastian_gsaenger@web.de",
    url="https://github.com/hein09/vipster",
    requires=["numpy", "OpenGL (>3.1.0)", "PyQt4"],
    scripts=["scripts/vipster"],
    packages=["vipster", "vipster.ioplugins",
              "vipster.gui", "vipster.gui.tools"],
    package_data={"vipster": ["default.json"],
                  "vipster.gui": ["opengl/*.vert", "opengl/*.frag",
                                  "opengl/bond_model", "opengl/sphere_model"]},
    data_files=[("share/doc/vipster", ["LICENSE", "README.md"]),
                ("share/pixmaps", ["vipster.png", "vipster-icon.png"]),
                ("share/applications", ["vipster.desktop"])],
    ext_modules=[mol_c, gui_c],
    license="BSD",
    classifiers=["Development Status :: 4 - Beta",
                 "Environment :: X11 Applications :: Qt",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: BSD License",
                 "Natural Language :: English",
                 "Programming Language :: Python :: Implementation :: CPython",
                 "Programming Language :: C",
                 "Topic :: Scientific/Engineering :: Chemistry",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Visualization"]
)
