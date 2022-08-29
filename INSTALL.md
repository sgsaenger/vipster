# Installation instructions

Vipster has four main components:

- *QtVipster*, a desktop-application
- *PyVipster*, python bindings to accomodate scripting
- *WebVipster*, a web-application
- *libvipster*, the library everything else is based on

## Precompiled releases

Binary releases of *QtVipster* are provided for Windows, Linux and (infrequently) macOS on [github](https://github.com/sgsaenger/vipster/releases).
Preliminary builds are provided as artifacts in the [CI pipelines](https://github.com/sgsaenger/vipster/actions) (github login required).

An exemplary implementation of *WebVipster* can be found [here](https://sgsaenger.github.io/vipster/emscripten).

An easy way to install *PyVipster* is via pip and PyPi: `pip install vipster`

### QtVipster on Windows

Windows releases are distributed as .7z archives, so please make sure to unpack it with [7-Zip](https://7-zip.org).

### QtVipster on Linux

For Arch Linux, PKBUILDs are available in the AUR as *vipster* (latest release) and *vipster-git* (following the latest developments).
This also contains the python bindings.

Other distributions may follow.

If your distribution does not offer a package, there is a portable .AppImage file you can download from the releases page. This can be made executable and used directly:

```sh
chmod +x Vipster-Linux-x86_64.AppImage
./Vipster-Linux-x86_64.AppImage
```

## Build from source

In order to compile Vipster, you need a working [CMake (>= 3.14)](https://cmake.org) installation.

Your C++ compiler should support C++17.
Vipster is tested against GCC (>=8) and Clang (>=4).
So far, MSVC is not supported, please use MinGW (>=9) if you are on Windows (see e.g. [MSYS2](https://www.msys2.org)).

### Dependencies

*QtVipster* requires Qt in version >=5.10.
For building from source, the Qt header files are required.
On Windows, they are included with the default installation.
On Linux, you may be required to install an additional package, e.g. `qtbase5-dev` or `qt6-base-dev` on Ubuntu.

Other dependencies are listed in `external/README.md`.
If they are not installed in your system,
they are included as git submodules and will be configured automatically during the build.
(This can be disabled by providing `-DVIPSTER_DOWNLOAD_DEPENDENCIES=OFF` to CMake, see below.)

### Build process (simple)

We need two directories:

- `$BUILD_DIR` is the directory that will contain your compiled files and can be chosen freely
- `$SOURCE_DIR` is the path to the Vipster source tree, e.g. where you unzipped the download or cloned the Git tree.

All commands shall be executed in `$BUILD_DIR`.
Vipster is built in multiple steps:

Step 1: Configure the build via CMake (see below for more options):

```sh
cmake $SOURCE_DIR
```

Step 2: Perform the compilation:

```sh
cmake --build .
```

Step 3: (Optional, not Windows) Install your compiled program:

```sh
cmake --build . -- install
```

### Build options (advanced)

You can provide options to CMake to further configure your installation.
All options are given in the form of `-D<OPTION>=<VALUE>`.
Step 1 can be repeated multiple times to apply modifications to the current configuration.
You can also use `cmake-gui` (graphical) or `ccmake` (terminal) to change options interactively.
After changing the configuration, execute Step 2 (and 3) again to update your compiled (installed) files.

List of common/vipster options, see [CMake documentation](https://cmake.org/cmake/help/latest/manual/cmake-variables.7.html) for more information:

- CMAKE_BUILD_TYPE: General compilation preset, set to "Release" for a regular optimized build or "Debug" for debug flags
- CMAKE_INSTALL_PREFIX: Set install path, defaults to `/usr/local`. Common setting for a per-user install: `$HOME/.local`
- VIPSTER_DESKTOP: Set to "ON" to build QtVipster, "OFF" otherwise (default if Qt is found)
- VIPSTER_LAMMPS: Set to "ON" to enable embedded LAMMPS widget (depends on QtVipster).
                  Will use a suitable installation if found, else will build it in-tree.
                  See the [LAMMPS-Manual](https://lammps.sandia.gov/doc/Manual.html) for more information and configuration options.
- VIPSTER_PYWIDGET: Set to "ON" to enable embedded Python widget (depends on QtVipster, defaults to "ON" if Python is found)
- VIPSTER_DOWNLOAD_DEPENDENCIES: Set to "OFF" to disable automatic git submodule initialization
- BUILD_TESTING: Set to "OFF" to disable testing ("ON" by default)
- CMAKE_PREFIX_PATH: Set to path for your Qt installation if not found automatically

### Debugging/Tests

If you intend to debug this software, you can enable debug information and disable optimizations by specifying `-D CMAKE_BUILD_TYPE=Debug`.

Unit tests are built by default (see `BUILD_TESTING`) and can be executed via `ctest`.

### PyVipster

To build the python library from source, execute the `setup.py` script in the root folder:

```sh
python setup.py install --user # will install into your home directory
python setup.py bdist_wheel # will create a .wheel file for distribution
```

#### For package maintainers

Installation via `setup.py` or `pip` will provide a statically linked library.
To get bindings that dynamically link to a system-wide installation,
use the CMake flag `-DVIPSTER_PYLIB=ON`.

### Web-Frontend

A feature-reduced browser frontend is available for browsers supporting WebGL2 and WebAssembly.
Use [emscripten](http://kripken.github.io/emscripten-site) for building:

```sh
cd $BUILD_DIR
emcmake cmake -D CMAKE_BUILD_TYPE=Release $SOURCE_DIR
emcmake cmake --build .
```

This prepares a .wasm file that contains the code, and a .js file that contains the Javascript interface.
To use this, one needs to embed it in a webpage and bind the exposed functions to HTML-events.
An example implementation can be found in `gh-pages/emscripten`.

#### Rebuild CSS

To rebuild the webpage's css, run

```sh
npm install
```

in `gh-pages/emscripten`. This will pull in dependencies and rebuild `styles.css`.

The following NPM scripts are available:

- `npm run build`:

    Compiles the SCSS code and executes the PostCSS Autoprefixer afterwards

- `npm run css:watch`:

    Starts the watch-mode: As soon as one of the SCSS files under `styles/` is changed, a CSS rebuild is triggered automatically.
