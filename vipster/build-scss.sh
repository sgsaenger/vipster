#!/usr/bin/env bash

SCRIPT=`realpath $0`
PWD=`dirname $SCRIPT`

sassc -t compressed $PWD/styles/styles.scss $PWD/styles/styles.css
npx postcss $PWD/styles/styles.css --use autoprefixer -r
