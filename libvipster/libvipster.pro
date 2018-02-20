QT       -= core gui

TARGET = vipster
TEMPLATE = lib

DEFINES += LIBVIPSTER_LIBRARY

CONFIG += c++14

SOURCES += \
    step.cpp \
    stepproper.cpp \
    stepformatter.cpp \
    molecule.cpp \
    config.cpp \
    iowrapper.cpp \
    ioplugins/xyz.cpp \
    ioplugins/pwinput.cpp \
    ioplugins/lmpinput.cpp \
    ioplugins/lmptrajec.cpp \
    atom.cpp \
    atomproper.cpp \
    atomref.cpp \
    ioplugins/pwoutput.cpp

HEADERS += \
    config.h \
    global.h \
    atom.h \
    bond.h \
    vec.h \
    kpoints.h \
    molecule.h \
    step.h \
    stepproper.h \
    stepformatter.h \
    iowrapper.h \
    ioplugin.h \
    ioplugins/xyz.h \
    ioplugins/pwinput.h \
    ioplugins/lmpinput.h \
    ioplugins/lmptrajec.h \
    atomproper.h \
    atomref.h \
    ioplugins/pwoutput.h

win32: CONFIG += staticlib
unix {
    target.path = /usr/lib
    INSTALLS += target
}
unix:!macx:CONFIG(debug, debug|release) {
    LIBS += -lgcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage -O0
    QMAKE_LFLAGS += -fprofile-arcs -ftest-coverage -O0
}

DISTFILES += \
    default.json
