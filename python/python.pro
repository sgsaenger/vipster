#-------------------------------------------------
#
# Project created by QtCreator 2017-01-10T00:29:28
#
#-------------------------------------------------

QT       -= core gui

TARGET = vipster
TEMPLATE = lib

DEFINES += LIBVIPSTER_PYTHON
CONFIG += c++14 plugin no_plugin_name_prefix

SOURCES += \
    vipster.cpp

HEADERS +=

QMAKE_CFLAGS += $$system(python-config --cflags)
QMAKE_CXXFLAGS += $$system(python-config --cflags)

unix {
    target.path = /usr/lib
    INSTALLS += target
}

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../libvipster/release/ -lvipster
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../libvipster/debug/ -lvipster
else:unix: LIBS += -L$$OUT_PWD/../libvipster/ -lvipster

INCLUDEPATH += $$PWD/../libvipster
DEPENDPATH += $$PWD/../libvipster
