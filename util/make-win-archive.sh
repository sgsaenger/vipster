#!/usr/bin/bash
# create a zip file from a pre-built tree

echo "creating zip-Archive"
TMPDIR=deploy
mkdir ${TMPDIR}

# Vipster itself
echo "Copying vipster"
cp gui/qt/vipster.exe ${TMPDIR}
cp vipster/vipster.dll ${TMPDIR}
cp vipster/vipster.lib ${TMPDIR}

# copy the local python installation
echo "Obtain python"
PY_LIB=$(grep "Python_LIBRARY_RELEASE" CMakeCache.txt | cut -d "=" -f 2)
cp $PY_LIB ${TMPDIR}
cp -r $pythonLocation/Lib $pythonLocation/DLLs ${TMPDIR}

# prepare Qt and other libraries
echo "Run windeployqt"
windeployqt --compiler-runtime --no-translations ${TMPDIR}/vipster.exe

echo "Create zip archive"
7z a Vipster-Win-x86_64.zip ${TMPDIR}/*
