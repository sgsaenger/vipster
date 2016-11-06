#-------------------------------------------------
#
# Project created by QtCreator 2016-07-27T14:00:20
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = vipster
TEMPLATE = app

CONFIG += c++11

SOURCES += main.cpp\
        mainwindow.cpp \
    glwidget.cpp \
    molwidget.cpp

HEADERS  += mainwindow.h \
    glwidget.h \
    atom_model.h \
    bond_model.h \
    molwidget.h

FORMS    += mainwindow.ui \
    molwidget.ui

DISTFILES += \
    resources/atom.frag \
    resources/atom.vert \
    resources/bond.frag \
    resources/bond.vert \
    resources/cell.frag \
    resources/cell.vert \
    resources/icon.png

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../libvipster/release/ -lvipster
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../libvipster/debug/ -lvipster
else:unix: LIBS += -L$$OUT_PWD/../libvipster/ -lvipster
unix:CONFIG(debug, debug|release) {
    LIBS += -lgcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage -O0
    QMAKE_LFLAGS += -fprofile-arcs -ftest-coverage -O0
}

INCLUDEPATH += $$PWD/../libvipster
DEPENDPATH += $$PWD/../libvipster

RESOURCES += \
    resources/vipster.qrc
