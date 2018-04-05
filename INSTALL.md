# Installation instructions

Vipster has three main frontends:
- QtVipster, a desktop-application
- PyVipster, python bindings to accomodate scripting
- WebVipster, a web-application

In order to compile any of these (except the web-frontend),
you need a suitable **qbs**-profile.
If you haven't configured your **qbs**, a straight-forward way is this:

```
qbs setup-toolchains --detect
qbs setup-qt --detect
```
Please pay attention to the output of these commands.
**qbs** will try to attach a toolchain to your qt installation, but if it fails,
you will have to execute this additional command:
```
qbs config profiles.$QT.baseProfile $TOOLCHAIN
```
where `$QT` and `$TOOLCHAIN` shall be the auto-generated names of the setup-X calls.


## Qt-Frontend (Win, Lin)

This frontend works on every pc with OpenGL3.3 capabilities (MacOS X postponed due to unknown problems).
To build the frontend, run:

```
cd $ARBITRARY_BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:$QT qbs.buildVariant:$VARIANT
```

where `$VIPSTER_SOURCE` shall be the directory you cloned the git/unpacked the archive in.

`$VARIANT` can be one of debug, release, or profile.
Debug will be the default, so please set it to release if you just want a working program.


## Python-bindings (Win, Lin, MacOS X)

These bindings should work on every platform that has python in its `$PATH`.
For Windows and MacOS, however, no installation is performed,
so you have to make sure that both the Vipster-library *and* the python-module are where they should be.
Maybe this will be fixed somewhen.
To build it, run:

```
cd $ARBITRARY_BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:$QT qbs.buildVariant:$VARIANT project.pythonBuild:true
```


## Web-Frontend

This frontend works in every browser with WebGL2 and WebAssembly support.
The needed qbs-profile is included in the project, so no need for setup.
It expects emcc(++) in your `$PATH`, so please make sure it is properly set.
To build the frontend, run:

```
cd $ARBITRARY_BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:emscripten qbs.buildVariant:release project.webBuild:true
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
