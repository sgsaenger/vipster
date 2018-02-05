#!/usr/bin/env bash
SCRIPT=`realpath $0`
PWD=`dirname $SCRIPT`

sassc -t compressed --sourcemap $PWD/styles/styles.scss $PWD/styles/styles.css
