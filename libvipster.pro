#-------------------------------------------------
#
# Project created by QtCreator 2016-08-03T17:38:47
#
#-------------------------------------------------

QT       -= core gui

TARGET = vipster
TEMPLATE = lib

DEFINES += LIBVIPSTER_LIBRARY

CONFIG += c++14

SOURCES += \
    config.cpp \
    molecule.cpp \
    step.cpp \
    definitions.cpp \
    iowrapper.cpp \
    ioplugins/xyz.cpp

HEADERS += \
    config.h \
    molecule.h \
    json.hpp \
    definitions.h \
    libvipster.h \
    step.h \
    iowrapper.h \
    ioplugins/xyz.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DISTFILES += \
    default.json
