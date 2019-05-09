# [![vipster](util/vipster.png)](https://sgsaenger.github.io/vipster) VIsual Periodic STructure EditoR

[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/2166/badge)](https://bestpractices.coreinfrastructure.org/projects/2166)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a276a159c93f47768c59dc264750f9f5)](https://www.codacy.com/app/sgsaenger/vipster?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=sgsaenger/vipster&amp;utm_campaign=Badge_Grade)

Master:
[![Build Status](https://travis-ci.com/sgsaenger/vipster.svg?branch=master)](https://travis-ci.com/sgsaenger/vipster)
[![Build Status (Windows)](https://ci.appveyor.com/api/projects/status/caoyp2efkyt6ly3x/branch/master?svg=true)](https://ci.appveyor.com/project/sgsaenger/vipster/branch/master)
[![codecov](https://codecov.io/gh/sgsaenger/vipster/branch/master/graph/badge.svg)](https://codecov.io/gh/sgsaenger/vipster)

Testing:
[![Build Status](https://travis-ci.com/sgsaenger/vipster.svg?branch=testing)](https://travis-ci.com/sgsaenger/vipster)
[![Build Status (Windows)](https://ci.appveyor.com/api/projects/status/caoyp2efkyt6ly3x/branch/testing?svg=true)](https://ci.appveyor.com/project/sgsaenger/vipster/branch/testing)
[![codecov](https://codecov.io/gh/sgsaenger/vipster/branch/testing/graph/badge.svg)](https://codecov.io/gh/sgsaenger/vipster)

Visualization and editing framework for atomistic simulations.

For more information, please visit the [Homepage](https://sgsaenger.github.io/vipster).

Binary releases are available [here](https://github.com/sgsaenger/vipster/releases),
installation instructions can be found [here](INSTALL.md).

Most importantly, [try it in your browser!](https://sgsaenger.github.io/vipster/emscripten/index.html)

<table align="center">
  <tr>
    <th colspan=3>
      <img src="util/vipster.png" height=16>
      Libvipster
    </th>
  </tr>
  <tr>
    <td colspan=3>C++17 based backbone: Powerful container-classes and I/O</td>
  </tr>
  <tr>
    <th>
      <img src="https://s3-eu-west-1.amazonaws.com/qt-files/logos/built-with-Qt_Horizontal_Small.png" alt="Qt GUI" height=18>
    </th>
    <th>
      <img src="https://raw.githubusercontent.com/emscripten-core/emscripten/master/media/switch_logo.png" alt="Emscripten port" height=60>
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
- [CLI11](https://github.com/CLIUtils/CLI11) (included)
- [tinyexpr](https://github.com/codeplea/tinyexpr) (included)
- [CMake >= 3.9](https://cmake.org)
- and a C++17-capable compiler (g++/mingw > 7 or clang > 4)
- optional:
    - [Qt5 >= 5.7](https://www.qt.io) (desktop application)
    - [emscripten](http://kripken.github.io/emscripten-site) (web interface)
    - [pybind11](https://github.com/pybind/pybind11) (script interface, included)
    - [Catch2](https://github.com/catchorg/Catch2) (testing, included)

## Supported file types:

| Format          | Reading  | Writing  |
|-----------------|----------|----------|
| xyz (augmented) | &#10004; | &#10004; |
| PWScf input     | &#10004; | &#10004; |
| PWScf output    | &#10004; |          |
| LAMMPS data     | &#10004; | &#10004; |
| LAMMPS dump     | &#10004; |          |
| CPMD input      | &#10004; | &#10004; |
| Gaussian cube   | &#10004; |          |
| XCrysden        | &#10004; |          |
| ORCA input      | &#10004; | &#10004; |
