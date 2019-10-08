#!/usr/bin/bash
# create a dmg file from a pre-built tree

mkdir -p vipster.app/Contents/Frameworks
cp -a vipster.framework vipster.app/Contents/Frameworks
# fix rpath
export VIPVER=$(ggrep "Vipster VERSION" $TRAVIS_BUILD_DIR/CMakeLists.txt | ggrep -o "[0-9.]*")
install_name_tool -change @rpath/vipster.framework/Versions/$VIPVER/vipster @executable_path/../Frameworks/vipster.framework/Versions/$VIPVER/vipster vipster.app/Contents/MacOS/vipster
# create .dmg file
/usr/local/opt/qt/bin/macdeployqt vipster.app -dmg
mv vipster.dmg Vipster-OSX.dmg
