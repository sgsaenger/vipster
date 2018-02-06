# Vipster Web Frontend

## Install

In current directory `vipster/vipster`:

```
npm install
```

```
VIPSTER_SOURCE=~/git/vipster
WEB_BUILD=~/foo
mkdir -p $WEB_BUILD && cd $_
qmake -spec $VIPSTER_SOURCE/emscriptenmkspec $VIPSTER_SOURCE/vipster.pro
make
```

(Adapt `VIPSTER_SOURCE` and `WEB_BUILD` to your needs.)

## Usage

The following NPM scripts are available:

* `npm run build`:

    Compiles the SCSS code and executes the PostCSS Autoprefixer afterwards

* `npm run css:watch`:

    Starts the watch-mode: As soon as one of the SCSS files under `styles/` is changed, a CSS rebuild is triggered automatically.

**Don't forget to run the `make` command in the `$WEB_BUILD` directory to copy the CSS and such to the build directory finally!**
