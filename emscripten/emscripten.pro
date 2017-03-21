TEMPLATE = app

CONFIG += c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = emcc
QMAKE_CXXFLAGS = --bind -s USE_WEBGL2=1
QMAKE_LINK = emcc
QMAKE_LFLAGS = --bind --embed-file ../libvipster/default.json@vipster.json

INCLUDEPATH += $$PWD/../libvipster $$PWD/../vipster /usr/lib/emscripten/system/include

TARGET = vipster.html

SOURCES += \
    vipster.cpp \
    $$PWD/../libvipster/config.cpp \
    $$PWD/../libvipster/molecule.cpp \
    $$PWD/../libvipster/step.cpp
