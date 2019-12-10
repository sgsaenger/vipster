#!/usr/bin/bash
# create a 7z file from a pre-built tree

echo "creating 7z-Archive"
mkdir Vipster
cp vipster.exe Vipster
cp libvipster.dll Vipster
cp $pythonLocation/python38.dll Vipster
7z a -tzip Vipster/python38.zip $pythonLocation/Lib/*
windeployqt --release --compiler-runtime Vipster/vipster.exe
7z a Vipster-Win-x86_64.7z Vipster
