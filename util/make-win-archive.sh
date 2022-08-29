#!/usr/bin/bash
# create a zip file from a pre-built tree

echo "creating zip-Archive"
mkdir Vipster

# Vipster itself
cp vipster.exe Vipster
cp libvipster.dll Vipster

# copy the local python installation
PY_LIB=$(grep "Python3_LIBRARY_RELEASE" CMakeCache.txt | cut -d "=" -f 2)
cp $PY_LIB Vipster
cp -r $pythonLocation/Lib $pythonLocation/DLLs Vipster

# prepare Qt and other libraries
windeployqt --compiler-runtime --no-translations Vipster/vipster.exe
cp $Qt5_Dir/bin/lib{gcc_s_seh-1,stdc++-6,winpthread-1}.dll Vipster

7z a Vipster-Win-x86_64.zip Vipster
