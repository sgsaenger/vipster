# [![logo](util/vipster.png)](https://sgsaenger.github.io/vipster) VIsual Periodic STructure EditoR

Master branch:
![Build status (master)](https://github.com/sgsaenger/vipster/workflows/Build/badge.svg?branch=master)

Testing branch:
![Build status (testing)](https://github.com/sgsaenger/vipster/workflows/Build/badge.svg?branch=testing)

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/sgsaenger/vipster)](https://github.com/sgsaenger/vipster/releases)
[![PyPI version](https://img.shields.io/pypi/v/vipster)](https://pypi.org/project/vipster)
[![Python versions](https://img.shields.io/pypi/pyversions/vipster)](https://pypi.org/project/vipster)

[![codecov](https://codecov.io/gh/sgsaenger/vipster/branch/master/graph/badge.svg)](https://codecov.io/gh/sgsaenger/vipster)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/2166/badge)](https://bestpractices.coreinfrastructure.org/projects/2166)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a276a159c93f47768c59dc264750f9f5)](https://app.codacy.com/manual/sgsaenger/vipster/dashboard)

[![DOI](https://zenodo.org/badge/21859848.svg)](https://zenodo.org/badge/latestdoi/21859848)
[![GPL-3.0 licensed](https://img.shields.io/github/license/sgsaenger/vipster)](https://www.gnu.org/licenses/gpl-3.0.html)

Fast and easy to use graphical editor for periodic atomistic simulations.
For more information, please visit the [Homepage](https://sgsaenger.github.io/vipster).
![Example screenshot](gh-pages/images/screenshot.png)

Binary releases are available [here](https://github.com/sgsaenger/vipster/releases),
installation instructions can be found [here](INSTALL.md).

Last but not least, a (sadly, feature-reduced) browser version can be tried [here](https://sgsaenger.github.io/vipster/emscripten/index.html)!

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
    </th>
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

- [JSON for Modern C++](https://github.com/nlohmann/json) 
- [CLI11](https://github.com/CLIUtils/CLI11) 
- [tinyexpr](https://github.com/codeplea/tinyexpr) 
- [CMake](https://cmake.org)
- [{fmt}](https://github.com/fmtlib/fmt)
- and a C++17-capable compiler (g++ > 8 or clang > 4)
- optional:
    - [Qt5](https://www.qt.io) (desktop application)
    - [emscripten](http://kripken.github.io/emscripten-site) (web interface)
    - [pybind11](https://github.com/pybind/pybind11) (script interface)
    - [Catch2](https://github.com/catchorg/Catch2) (testing)
    - [LAMMPS](https://lammps.sandia.gov) (interactive forcefield calculations)

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
| VASP Poscar     | &#10004; | &#10004; |
