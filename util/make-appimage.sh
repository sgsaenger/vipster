#!/usr/bin/env bash
# create a portable AppImage from a pre-built tree

echo "creating .AppImage"
# install vipster into AppDir
make DESTDIR=AppDir install

# extract build configuration from CMakeCache.txt
SOURCE_DIR=$(grep Vipster_SOURCE_DIR CMakeCache.txt | cut -d "=" -f 2)
PY_BIN=$(grep "Python3_EXECUTABLE" CMakeCache.txt | cut -d "=" -f 2)
PY_PREFIX=$(${PY_BIN} -c "import sys; print(sys.prefix)")

# install hook for dynamic selection of python runtime
mkdir -p AppDir/apprun-hooks
sed "s|X_PYROOT_X|${PY_PREFIX}|" "${SOURCE_DIR}/util/select_python_hook.sh" | sed "s|X_PYBIN_X|${PY_BIN##*/}|" > AppDir/apprun-hooks/select_python_hook.sh

# use linuxdeploy to manage dependencies
wget -q https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
wget -q https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/continuous/linuxdeploy-plugin-qt-x86_64.AppImage
wget -q https://github.com/linuxdeploy/linuxdeploy-plugin-checkrt/releases/download/continuous/linuxdeploy-plugin-checkrt-x86_64.sh
wget -q https://github.com/linuxdeploy/linuxdeploy-plugin-appimage/releases/download/continuous/linuxdeploy-plugin-appimage-x86_64.AppImage
chmod +x linuxdeploy*
./linuxdeploy-x86_64.AppImage --appdir AppDir -p qt -p checkrt

# bundle the python installation used for building
PY_LIB_DIR=$(${PY_BIN} -c "import sysconfig as s; print(s.get_path('stdlib'))")
mkdir -p AppDir/$PY_LIB_DIR
cp -r $PY_LIB_DIR/. AppDir/$PY_LIB_DIR

# move libpython so we can use system's version if compatible
mkdir -p AppDir/usr/optional/python
mv AppDir/usr/lib/libpython* AppDir/usr/optional/python

# finish by creating AppImage
ARCH=x86_64 OUTPUT=Vipster-Linux-x86_64.AppImage ./linuxdeploy-plugin-appimage-x86_64.AppImage --appdir AppDir
