#-------------------------------------------------
#
# Project created by QtCreator 2016-07-27T14:00:20
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = vipster
TEMPLATE = app

CONFIG += c++14

SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    atom_model.h \
    bond_model.h

FORMS    += mainwindow.ui

DISTFILES += \
    vipster-icon.png \
    atom.frag \
    atom.vert \
    bond.frag \
    bond.vert \
    cell.frag \
    cell.vert

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../libvipster/release/ -lvipster
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../libvipster/debug/ -lvipster
else:unix: LIBS += -L$$OUT_PWD/../libvipster/ -lvipster

INCLUDEPATH += $$PWD/../libvipster
DEPENDPATH += $$PWD/../libvipster

RESOURCES += \
    vipster.qrc
