#!/usr/bin/bash

# build release-version
cd $TRAVIS_BUILD_DIR
mkdir release
cd release
cmake -DDESKTOP=ON -DPYTHON=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr $TRAVIS_BUILD_DIR
make DESTDIR=AppDir install -j2
# copy python standard-library and add libpython to LD_LIBRARY_PATH
export PY_LIB_DIR=$(python -c "from distutils import sysconfig as s; print(s.get_python_lib(standard_lib=True))")
cp -r $PY_LIB_DIR/. AppDir/$PY_LIB_DIR
export LD_LIBRARY_PATH=$(python -c "from distutils import sysconfig as s; print(s.get_config_var('LIBDIR'))")
# fill AppDir with dependencies
wget -c "https://github.com/probonopd/linuxdeployqt/releases/download/6/linuxdeployqt-6-x86_64.AppImage" -O linuxdeployqt
chmod +x linuxdeployqt
./linuxdeployqt AppDir/usr/share/applications/vipster.desktop -bundle-non-qt-libs;
# move libpython so we can use possibly use system's version
mkdir -p AppDir/usr/optional/python
mv AppDir/usr/libpython* AppDir/usr/optional/python
# bundle libstdc++ to be compatible with older linuxes, see https://github.com/darealshinji/AppImageKit-checkrt
mkdir -p AppDir/usr/optional/libstdc++
wget -c https://github.com/darealshinji/AppImageKit-checkrt/releases/download/continuous/exec-x86_64.so -O AppDir/usr/optional/exec.so
cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 AppDir/usr/optional/libstdc++/
# use custom AppRun so we get the correct working directory
rm AppDir/AppRun
sed "s|X_PYROOT_X|$(pyenv prefix)|" $TRAVIS_BUILD_DIR/util/AppRun > AppDir/AppRun
chmod a+x AppDir/AppRun
# create AppImage
wget -c https://github.com/AppImage/AppImageKit/releases/download/11/appimagetool-x86_64.AppImage -O appimagetool
chmod +x appimagetool
./appimagetool -g AppDir $TRAVIS_BUILD_DIR/Vipster-Linux-x86_64.AppImage
