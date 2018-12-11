# Installation instructions

Vipster has four main components:
- *QtVipster*, a desktop-application
- *PyVipster*, python bindings to accomodate scripting
- *WebVipster*, a web-application
- *libvipster*, the library everything else is based on

In order to compile any of these, you need a working **cmake** installation.

## Qt-Frontend

This frontend should work on every pc with OpenGL3.3 capabilities.
Please make sure that you have a valid Qt installation.
To build the frontend, run:

```
cd $VIPSTER_SOURCE
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release -D DESKTOP=YES ..
cmake --build .
```

If cmake does not find your Qt installation (e.g. `qmake` not in your `$PATH`) or you want to specify a certain version,
add `-D CMAKE_PREFIX_PATH=$QT_ROOT` to the first cmake call.

## Python-bindings

These bindings should work on every platform that has python in its `$PATH`.
Please refer to the documentation of [pybind11](https://github.com/pybind/pybind11)

```
cd $VIPSTER_SOURCE
git submodule update --init --recursive
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release -D PYTHON=YES ..
cmake --build .
```

The python-bindings and Qt-frontend can be enabled at once by specifying both `-D PYTHON=YES -D DESKTOP=YES`

## Debugging

If you intend to debug this software, you can enable debug information and disable optimizations by specifying `-D CMAKE_BUILD_TYPE=Debug`.
Additionally, you can enable compilation of unit tests by adding `-D TESTS=YES`.

## Web-Frontend

This frontend works in every browser with WebGL2 and WebAssembly support.
It expects a working **emscripten** installation.
To compile, use:
```
cd $VIPSTER_SOURCE
mkdir build
cd build
emcmake cmake -D CMAKE_BUILD_TYPE=Release -D WEB=YES -D DESKTOP=NO ..
emcmake cmake --build .
```

This prepares a .wasm file that contains the code, and a .js file that contains the Javascript interface.
To use this, one needs to embed it in a webpage and bind the exposed functions to HTML-events.
An example implementation can be found in `gh-pages/emscripten`.

### Rebuild CSS (optional)

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
