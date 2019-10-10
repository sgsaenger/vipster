# Installation instructions

Vipster has four main components:
- *QtVipster*, a desktop-application
- *PyVipster*, python bindings to accomodate scripting
- *WebVipster*, a web-application
- *libvipster*, the library everything else is based on

## Precompiled releases

Binary releases of *QtVipster* are provided for Windows, Linux and macOS on [github](https://github.com/sgsaenger/vipster/releases).
An exemplary implementation of *WebVipster* can be found [here](https://sgsaenger.github.io/vipster/emscripten).
An easy way to install *PyVipster* is via pip and PyPi: `pip install vipster`

### QtVipster on Windows

Windows releases are distributed as .7z archives, so please make sure to unpack it with [7-Zip](https://7-zip.org).

### QtVipster on Linux

For Arch Linux, PKBUILDs are available in the AUR as *vipster* (latest release) and *vipster-git* (following the latest developments).
This also contains the python bindings.

Other distributions may follow.

If your distribution does not offer a package, there is a portable .AppImage file you can download from the releases page. This can be made executable and used directly:
```
chmod +x Vipster-Linux-x86_64.AppImage
./Vipster-Linux-x86_64.AppImage
```

### QtVipster on macOS

The macOS release is distributed as a regular .dmg file that you can install as usual.

## Build from source

In order to compile Vipster, you need a working [cmake](https://cmake.org) installation.
Your C++ compiler should support C++17.
Vipster is tested against GCC and Clang.
So far, MSVC is not supported, please use MinGW if you are on Windows.

### Qt-Frontend

This frontend should work on every pc with OpenGL3.3 capabilities.
Please make sure that you have a valid Qt installation.
To build the frontend, run:

```
cd $VIPSTER_SOURCE
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release -D DESKTOP=ON ..
cmake --build .
```

One may also want to add `make install` to deploy the files.
The target directory for this is set by adding e.g. `-D CMAKE_INSTALL_PREFIX=/usr`.
On linux, the default target directory will probably default to `$HOME/.local`, which may suffice for a per-user install.

If cmake does not find your Qt installation (e.g. `qmake` not in your `$PATH`) or you want to specify a certain version,
add `-D CMAKE_PREFIX_PATH=$QT_ROOT` to the first cmake call.

#### Built-in Python-widget

By adding `-D PYSHELL=ON` to the build of the Qt-frontend, a python terminal with embedded vipster bindings is enabled.

### Standalone Python-bindings

These bindings should work on every platform that has python in its `$PATH`.

If you are on Linux and your distribution offers a Vipster-package, it should include the bindings.
Please also refer to the documentation of [pybind11](https://github.com/pybind/pybind11).
To build them from source, you can either do `python setup.py install`, which will install like any other python package, or go the manual way (and share the build environment/settings):

```
cd $VIPSTER_SOURCE
git submodule update --init --recursive
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release -D PYBIND=ON ..
cmake --build .
```
**NOTE**: By default, Python-Egg informations will be generated.
This can be disabled by specifying `-DEGG_INFO=OFF`.

### Debugging/Tests

If you intend to debug this software, you can enable debug information and disable optimizations by specifying `-D CMAKE_BUILD_TYPE=Debug`.
Additionally, you can enable compilation of unit tests by adding `-D TESTS=ON`.

### Web-Frontend

This frontend works in every browser with WebGL2 and WebAssembly support.
It expects a working [emscripten](http://kripken.github.io/emscripten-site) installation.
To compile, use:
```
cd $VIPSTER_SOURCE
mkdir build
cd build
emcmake cmake -D CMAKE_BUILD_TYPE=Release -D WEB=ON -D DESKTOP=OFF ..
emcmake cmake --build .
```

This prepares a .wasm file that contains the code, and a .js file that contains the Javascript interface.
To use this, one needs to embed it in a webpage and bind the exposed functions to HTML-events.
An example implementation can be found in `gh-pages/emscripten`.

#### Rebuild CSS (optional)

To rebuild the css, run
```
npm install
```
in `gh-pages/emscripten`. This will pull in dependencies and rebuild `styles.css`.

The following NPM scripts are available:

* `npm run build`:

    Compiles the SCSS code and executes the PostCSS Autoprefixer afterwards

* `npm run css:watch`:

    Starts the watch-mode: As soon as one of the SCSS files under `styles/` is changed, a CSS rebuild is triggered automatically.
