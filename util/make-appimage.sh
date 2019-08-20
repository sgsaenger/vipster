#!/usr/bin/bash

# build release-version
cd $TRAVIS_BUILD_DIR
mkdir release
cd release
cmake -D DESKTOP=YES -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=/usr $TRAVIS_BUILD_DIR
make DESTDIR=AppDir install -j2
# fill AppDir with dependencies
wget -c "https://github.com/probonopd/linuxdeployqt/releases/download/6/linuxdeployqt-6-x86_64.AppImage" -O linuxdeployqt
chmod +x linuxdeployqt
./linuxdeployqt AppDir/usr/share/applications/vipster.desktop -bundle-non-qt-libs;
# bundle libstdc++ to be compatible with older linuxes, see https://github.com/darealshinji/AppImageKit-checkrt
mkdir -p AppDir/usr/optional/libstdc++
wget -c https://github.com/darealshinji/AppImageKit-checkrt/releases/download/continuous/exec-x86_64.so -O AppDir/usr/optional/exec.so
cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 AppDir/usr/optional/libstdc++/
# use custom AppRun so we get the correct working directory
rm AppDir/AppRun
cp $TRAVIS_BUILD_DIR/util/AppRun AppDir
chmod a+x AppDir/AppRun
# create AppImage
wget -c https://github.com/AppImage/AppImageKit/releases/download/11/appimagetool-x86_64.AppImage -O appimagetool
chmod +x appimagetool
./appimagetool -g AppDir $TRAVIS_BUILD_DIR/Vipster-Linux-x86_64.AppImage
