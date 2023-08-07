#!/usr/bin/bash
# create a dmg file from a pre-built tree

echo "creating .dmg"
mkdir -p gui/qt/vipster.app/Contents/Frameworks
cp -a vipster/vipster.framework gui/qt/vipster.app/Contents/Frameworks
# fix rpath
export VIPVER=$(grep "CMAKE_PROJECT_VERSION:" ../CMakeLists.txt | cut -d "=" -f 2)
install_name_tool -change @rpath/vipster.framework/Versions/$VIPVER/vipster @executable_path/../Frameworks/vipster.framework/Versions/$VIPVER/vipster gui/qt/vipster.app/Contents/MacOS/vipster
# create .dmg file
macdeployqt gui/qt/vipster.app -dmg
mv vipster.dmg Vipster-macOS-x86_64.dmg
