# [![vipster](vipster-icon.png)](https://hein09.github.io/vipster) VIsual Periodic STructure EditoR

Master:
[![Build Status](https://travis-ci.org/hein09/vipster.svg?branch=master)](https://travis-ci.org/hein09/vipster)
[![Build status](https://ci.appveyor.com/api/projects/status/caoyp2efkyt6ly3x/branch/master?svg=true)](https://ci.appveyor.com/project/hein09/vipster/branch/master)
[![codecov](https://codecov.io/gh/hein09/vipster/branch/master/graph/badge.svg)](https://codecov.io/gh/hein09/vipster)

Testing:
[![Build Status](https://travis-ci.org/hein09/vipster.svg?branch=testing)](https://travis-ci.org/hein09/vipster)
[![Build status](https://ci.appveyor.com/api/projects/status/caoyp2efkyt6ly3x/branch/testing?svg=true)](https://ci.appveyor.com/project/hein09/vipster/branch/testing)
[![codecov](https://codecov.io/gh/hein09/vipster/branch/testing/graph/badge.svg)](https://codecov.io/gh/hein09/vipster)

Visualization and editing framework for atomistic simulations.

For more information, please visit the [Homepage](https://hein09.github.io/vipster),
installation instructions can be found [here](INSTALL.md).

Most importantly, [try it in your browser!](https://hein09.github.io/vipster/emscripten/index.html)

<table align="center">
  <tr>
    <th colspan=3>
      <img src="vipster-icon.png" height=16>
      Libvipster
    </th>
  </tr>
  <tr>
    <td colspan=3>C++14 based backbone: Powerful container-classes and I/O</td>
  </tr>
  <tr>
    <th>
      <img src="https://s3-eu-west-1.amazonaws.com/qt-files/logos/built-with-Qt_Horizontal_Small.png" alt="Qt GUI" height=18>
    </th>
    <th>
      <img src="https://github.com/kripken/emscripten/blob/master/media/switch_logo.png" alt="Emscripten port" height=60>
    <th>
      <img src="https://www.python.org/static/community_logos/python-logo-master-v3-TM.png" alt="Python bindings" height=36>
    </th>
  </tr>
  <tr>
    <td>Fast and native GUI with OGL3.3 based rendering</td>
    <td>Portable browser-based GUI, shared render-code</td>
    <td>Scripting interface for batch-processing</td>
  </tr>
</table>

## Dependencies:

- [JSON for Modern C++ >= 3.0](https://github.com/nlohmann/json) (included)
- [Qt5 >= 5.7](https://www.qt.io) including QBS >= 1.11
- and a C++14-capable compiler (g++/mingw > 5 or clang > 3.4)
- optional:
    - [Catch2](https://github.com/catchorg/Catch2) (testing, included)
    - [pybind11](https://github.com/pybind/pybind11) (script-interface, included)
    - [emscripten](http://kripken.github.io/emscripten-site) (web-interface)

## Supported file types:

- standard xyz
- PWscf input
- PWScf output
- Lammps data
- Lammps dump
- ~~Empire-xyz~~
- ~~Gaussian cube~~
- ~~AIMALL output~~
- ~~CPMD input~~
- ~~Mol2~~
- ~~XSF/AXSF~~
- ~~Turbomole~~

