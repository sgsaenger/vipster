# Vipster Frontend

## Qt-Frontend (Win, Lin)

This frontend works on every pc with OpenGL3.3 capabilities (MacOS X postponed due to unknown problems).
To build the frontend, run:

```
cd $WHEREVER
mkdir build-vipster && cd build-vipster
qmake $VIPSTER_SOURCE/vipster.pro
make
```
where `VIPSTER_SOURCE` shall be the directory you cloned the git/unpacked the archive in.

## Web-Frontend

This frontend works in every browser with WebGL2 and WebAssembly support.
To build the frontend, run:

```
cd $WHEREVER
mkdir build-webvipster && cd build-webvipster
qmake -spec $VIPSTER_SOURCE/emscriptenmkspec $VIPSTER_SOURCE/vipster.pro
make
```
where `VIPSTER_SOURCE` shall be the directory you cloned the git/unpacked the archive in.

### Rebuild CSS (optional)

To rebuild the css, run
```
npm install
```
in this directory. This will pull in dependencies and rebuild `styles/styles.css`.

The following NPM scripts are available:

* `npm run build`:

    Compiles the SCSS code and executes the PostCSS Autoprefixer afterwards

* `npm run css:watch`:

    Starts the watch-mode: As soon as one of the SCSS files under `styles/` is changed, a CSS rebuild is triggered automatically.

**Don't forget to run the `make` command in the `$WEB_BUILD` directory after every CSS change to copy the updates files to the build directory!**
