#-------------------------------------------------
#
# Project created by QtCreator 2017-01-10T00:29:28
#
#-------------------------------------------------

QT       -= core gui

TARGET = vipster
TEMPLATE = lib

DEFINES += LIBVIPSTER_PYTHON

SOURCES += \
    vipster.cpp

HEADERS +=

unix {
    target.path = /usr/lib
    INSTALLS += target
}
