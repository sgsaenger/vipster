TEMPLATE = app
CONFIG += c++14

INCLUDEPATH += $$PWD/../libvipster
DEPENDPATH += $$PWD/../libvipster

wasm {
    QT -= gui
    TARGET = vipster.js
    CONFIG -= app_bundle
    CONFIG -= qt

    SOURCES += guiwrapper.cpp \
        mainweb.cpp

    HEADERS += guiwrapper.h

    QMAKE_LFLAGS += \
        --embed-file $$PWD/../libvipster/default.json@vipster.json \
        --embed-file $$PWD/../vipster/resources/atom.frag@atom.frag \
        --embed-file $$PWD/../vipster/resources/atom.vert@atom.vert \
        --embed-file $$PWD/../vipster/resources/bond.frag@bond.frag \
        --embed-file $$PWD/../vipster/resources/bond.vert@bond.vert \
        --embed-file $$PWD/../vipster/resources/cell.frag@cell.frag \
        --embed-file $$PWD/../vipster/resources/cell.vert@cell.vert

} else {
    QT += core gui widgets
    TARGET = vipster

    SOURCES += main.cpp\
        mainwindow.cpp \
        glwidget.cpp \
        guiwrapper.cpp\
        molwidget.cpp

    HEADERS += mainwindow.h \
        glwidget.h \
        guiwrapper.h \
        atom_model.h \
        bond_model.h \
        molwidget.h

    FORMS += mainwindow.ui \
        molwidget.ui

    RESOURCES += \
        resources/vipster.qrc
}
win32: CONFIG += static no_smart_library_merge
win32:CONFIG(release, debug|release) {
    LIBS += -L$$OUT_PWD/../libvipster/release/ -lvipster
    RC_ICONS = resources/vipster.ico
} else:win32:CONFIG(debug, debug|release) {
    LIBS += -L$$OUT_PWD/../libvipster/debug/ -lvipster
} else:unix|wasm {
    LIBS += -L$$OUT_PWD/../libvipster/ -lvipster
}

unix:CONFIG(debug, debug|release) {
    LIBS += -lgcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage -O0
    QMAKE_LFLAGS += -fprofile-arcs -ftest-coverage -O0
}
