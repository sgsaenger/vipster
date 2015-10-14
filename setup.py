#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy as np

bonds = Extension(name="bonds",sources=["C/set_bonds_c.c"], include_dirs = [np.get_include()],extra_compile_args=['-std=c99'])

setup(
        name="pwtoolbox",
        version="0.5",
        description ="GUI/CLI interface for periodic structures",
        author ="Sebastian Gsänger",
        author_email="sebastian_gsaenger@web.de",
        url="https://github.com/hein09/pwtoolbox",
        requires=["numpy","OpenGL (>3.1.0)","PyQt4"],
        scripts=['ptb'],
        py_modules=[
            'ptb_mol','molecule',
            'collapsiblewidget','conftab','moltab','pwtab','viewport','ptb_gui'],
        packages=['tools','ftypeplugins'],
        ext_modules=[bonds]
)
