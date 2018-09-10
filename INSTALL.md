# Installation instructions

Vipster has three main frontends:
- QtVipster, a desktop-application
- PyVipster, python bindings to accomodate scripting
- WebVipster, a web-application

In order to compile any of these (except the web-frontend),
you need a suitable **qbs**-profile.
If you have not configured your **qbs**, a straight-forward way is this:

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


## Qt-Frontend

This frontend works on every pc with OpenGL3.3 capabilities (MacOS X postponed due to unknown problems).
To build the frontend, run:

```
cd $ARBITRARY_BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:$QT qbs.buildVariant:$VARIANT
```

where `$VIPSTER_SOURCE` shall be the directory you cloned the git/unpacked the archive in.

`$VARIANT` can be one of debug, release, or profile.
Debug will be the default, so please set it to release if you just want a working program.

### For Linux (and other Unices):
The default configuration expects to be installed under /usr/{bin,include,lib,share}.
If you would like to specify another install prefix, please add `project.prefix:"/target/prefix"` to the qbs line.
When installing in your home-directory (or for creating an AppImage or similar), set `project.relpath` to `true`.

Example for a user-specific install:
```
qbs build -f $VIPSTER_SOURCE profile:$QT qbs.buildVariant:release project.prefix:"/home/NAME/.local"
```

### Windows, MacOS X:

Adding `project.winInstall:true` or `project.macInstall:true` will create the respective installers.


## Python-bindings

These bindings should work on every platform that has python in its `$PATH`.
For Windows and MacOS, however, no installation is performed,
so you have to make sure that both the Vipster-library *and* the python-module are where they should be.
Maybe this will be fixed somewhen.
To build it, add `project.pythonBuild:true` to the qbs line.
If you would like to select a different python executable, add `project.pythonName:"/path/to/python"` to the qbs line.

```
cd $ARBITRARY_BUILD_DIR
qbs build -f $VIPSTER_SOURCE profile:$QT qbs.buildVariant:release project.pythonBuild:true project.pythonName:"python3"
```


## Web-Frontend

This frontend works in every browser with WebGL2 and WebAssembly support.
The needed qbs-profile (`emscripten`) is included in the project, so no need for setup.
It expects emcc(++) in your `$PATH`, so please make sure everything is setup correctly.
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

**Do not forget to rebuild the web-frontend after making some changes to the CSS!**
