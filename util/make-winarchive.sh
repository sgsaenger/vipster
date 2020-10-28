#!/usr/bin/bash
# create a 7z file from a pre-built tree

echo "creating 7z-Archive"
mkdir Vipster
cp vipster.exe Vipster
cp libvipster.dll Vipster
cp $pythonLocation/python39.dll Vipster
7z a -tzip Vipster/python39.zip $pythonLocation/Lib/*
windeployqt --release --compiler-runtime --no-translations Vipster/vipster.exe
cp $QTDIR/bin/lib{gcc_s_seh-1,stdc++-6,winpthread-1}.dll Vipster
7z a Vipster-Win-x86_64.7z Vipster
