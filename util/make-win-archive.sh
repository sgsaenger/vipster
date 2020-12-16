#!/usr/bin/bash
# create a zip file from a pre-built tree

echo "creating zip-Archive"
cd build
mkdir Vipster
cp vipster.exe Vipster
cp libvipster.dll Vipster
cp $pythonLocation/python38.dll Vipster
cp -r $pythonLocation/Lib $pythonLocation/DLLs Vipster
windeployqt --compiler-runtime --no-translations Vipster/vipster.exe
cp $Qt5_Dir/bin/lib{gcc_s_seh-1,stdc++-6,winpthread-1}.dll Vipster
7z a ../Vipster-Win-x86_64.zip Vipster
