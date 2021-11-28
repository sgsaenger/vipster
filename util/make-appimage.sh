#!/usr/bin/env bash
# create a portable AppImage from a pre-built tree

echo "creating .AppImage"
# install vipster into AppDir
make DESTDIR=AppDir install

# get source dir and correct python from CMakeCache.txt
export SOURCE_DIR=$(grep Vipster_SOURCE_DIR CMakeCache.txt | cut -d "=" -f 2)
export PY_BIN=$(grep "Python3_EXECUTABLE" CMakeCache.txt | cut -d "=" -f 2)
export PY_PREFIX=$(${PY_BIN} -c "from distutils import sysconfig as s; print(s.PREFIX)")

# add libpython to LD_LIBRARY_PATH so linuxdeployqt finds it
export LD_LIBRARY_PATH=$(${PY_BIN} -c "from distutils import sysconfig as s; print(s.get_config_var('LIBDIR'))"):$LD_LIBRARY_PATH

# bundle required libraries (Qt, Python)
wget -q -c "https://github.com/probonopd/linuxdeployqt/releases/download/6/linuxdeployqt-6-x86_64.AppImage" -O linuxdeployqt
chmod +x linuxdeployqt
./linuxdeployqt AppDir/usr/share/applications/com.github.sgsaenger.vipster.desktop -bundle-non-qt-libs -unsupported-allow-new-glibc

# copy python standard-library
export PY_LIB_DIR=$(${PY_BIN} -c "from distutils import sysconfig as s; print(s.get_python_lib(standard_lib=True))")
mkdir -p AppDir/$PY_LIB_DIR
cp -r $PY_LIB_DIR/. AppDir/$PY_LIB_DIR

# move libpython so we can use system's version if compatible
mkdir -p AppDir/usr/optional/python
mv AppDir/usr/lib/libpython* AppDir/usr/optional/python

# bundle libstdc++ to be compatible with older linuxes, see https://github.com/darealshinji/AppImageKit-checkrt
mkdir -p AppDir/usr/optional/libstdc++
wget -q -c https://github.com/darealshinji/AppImageKit-checkrt/releases/download/continuous/exec-x86_64.so -O AppDir/usr/optional/exec.so
cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 AppDir/usr/optional/libstdc++/

# use custom AppRun so we get the correct working directory and fallback libraries
rm AppDir/AppRun
sed "s|X_PYROOT_X|${PY_PREFIX}|" "${SOURCE_DIR}/util/AppRun" | sed "s|X_PYBIN_X|${PY_BIN##*/}|" > AppDir/AppRun
chmod a+x AppDir/AppRun

# bundle AppDir as AppImage
wget -q -c https://github.com/AppImage/AppImageKit/releases/download/11/appimagetool-x86_64.AppImage -O appimagetool
chmod +x appimagetool
./appimagetool -g AppDir Vipster-Linux-x86_64.AppImage
