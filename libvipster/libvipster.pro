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
    ioplugins/xyz.h \
    ioplugin.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
unix:CONFIG(debug, debug|release) {
    LIBS += -lgcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage -O0
    QMAKE_LFLAGS += -fprofile-arcs -ftest-coverage -O0
}

DISTFILES += \
    default.json
