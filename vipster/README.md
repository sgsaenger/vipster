# Vipster Frontend

## Qt-Frontend (Win, Lin)

This frontend works on every pc with OpenGL3.3 capabilities (MacOS X postponed due to unknown problems).
To build the frontend, run:

```
cd $BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:$YOUR_QT_PROFILE qbs.buildVariant:$VARIANT
```

where `VIPSTER_SOURCE` shall be the directory you cloned the git/unpacked the archive in.

`YOUR_QT_PROFILE` is the name of a qbs-profile you have to setup manually beforehand.
A simple working configuration can be yielded like so:
```
qbs setup-toolchains --detect
qbs setup-qt --detect
```
Please refer to the qbs-documentation for further information.

`VARIANT` can be one of debug, release, or profile, where you most likely want to choose release.

## Web-Frontend

This frontend works in every browser with WebGL2 and WebAssembly support.
To build the frontend, run:

```
cd $BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:emscripten qbs.buildVariant:release
```

where `VIPSTER_SOURCE` shall be the directory you cloned the git/unpacked the archive in.

### Rebuild CSS (optional)

To rebuild the css, run
```
npm install
```
in `vipster/web/page`. This will pull in dependencies and rebuild `styles.css`.

The following NPM scripts are available:

* `npm run build`:

    Compiles the SCSS code and executes the PostCSS Autoprefixer afterwards

* `npm run css:watch`:

    Starts the watch-mode: As soon as one of the SCSS files under `styles/` is changed, a CSS rebuild is triggered automatically.

**Don't forget to rebuild the web-frontend after making some changes to the CSS!**
