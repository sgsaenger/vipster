QT       -= core gui

TARGET = vipster
TEMPLATE = lib

DEFINES += LIBVIPSTER_LIBRARY

CONFIG += c++14

SOURCES += \
    step.cpp \
    molecule.cpp \
    config.cpp \
    iowrapper.cpp \
    ioplugins/xyz.cpp \
    ioplugins/pwinput.cpp \
    ioplugins/lmpinput.cpp \
    ioplugins/lmptrajec.cpp

HEADERS += \
    config.h \
    molecule.h \
    step.h \
    iowrapper.h \
    ioplugin.h \
    global.h \
    atom.h \
    bond.h \
    vec.h \
    kpoints.h \
    ioplugins/xyz.h \
    ioplugins/pwinput.h \
    ioplugins/lmpinput.h \
    ioplugins/lmptrajec.h

win32: CONFIG += staticlib
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
