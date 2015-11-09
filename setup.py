#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy as np

mol_c = Extension(name="ptb.mol_c",sources=["vipster/mol_c.c"], include_dirs = [np.get_include()],extra_compile_args=['-std=c99'])
gui_c = Extension(name="ptb.gui.gui_c",sources=["vipster/gui/gui_c.c"], include_dirs = [np.get_include()],extra_compile_args=['-std=c99'])

setup(
        name="pwtoolbox",
        version="0.8.0",
        description ="GUI for periodic structures",
        author ="Sebastian Gsänger",
        author_email="sebastian_gsaenger@web.de",
        url="https://github.com/hein09/pwtoolbox",
        requires=["numpy","OpenGL (>3.1.0)","PyQt4"],
        scripts=["ptb"],
        package_dir={"ptb":"vipster"},
        packages=["ptb","ptb.ftypeplugins","ptb.gui","ptb.gui.tools"],
        package_data={"ptb":["default.json"],"ptb.gui":["opengl/*.vert","opengl/*.frag","opengl/bond_model","opengl/sphere_model"]},
        data_files=[("share/doc/pwtoolbox",["LICENSE","README.md"])],
        ext_modules=[mol_c,gui_c],
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
