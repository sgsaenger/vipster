TEMPLATE = app

CONFIG += c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = emcc
QMAKE_CXXFLAGS = --bind
QMAKE_LINK = emcc
QMAKE_LFLAGS = --bind

INCLUDEPATH += $$PWD/../libvipster /usr/lib/emscripten/system/include

TARGET = vipster.html

SOURCES += \
    vipster.cpp
