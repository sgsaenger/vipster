TEMPLATE = app

CONFIG += c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = emcc
QMAKE_CXXFLAGS = --bind -s USE_WEBGL2=1
QMAKE_LINK = emcc
QMAKE_LFLAGS = --bind --embed-file $$PWD/../libvipster/default.json@vipster.json \
    --embed-file $$PWD/../vipster/resources/atom.frag@atom.frag \
    --embed-file $$PWD/../vipster/resources/atom.vert@atom.vert \
    --embed-file $$PWD/../vipster/resources/cell.frag@cell.frag \
    --embed-file $$PWD/../vipster/resources/cell.vert@cell.vert \
    -s USE_WEBGL2=1

INCLUDEPATH += $$PWD/../libvipster $$PWD/../vipster /usr/lib/emscripten/system/include

TARGET = vipster.js

SOURCES += \
    vipster.cpp \
    $$PWD/../libvipster/config.cpp \
    $$PWD/../libvipster/molecule.cpp \
    $$PWD/../libvipster/step.cpp \
    $$PWD/../libvipster/iowrapper.cpp \
    $$PWD/../libvipster/ioplugins/xyz.cpp \
    $$PWD/../libvipster/ioplugins/pwinput.cpp \
    $$PWD/../vipster/guiwrapper.cpp
