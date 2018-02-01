TEMPLATE = app
CONFIG += c++14

INCLUDEPATH += $$PWD/../libvipster
DEPENDPATH += $$PWD/../libvipster
ICON = resources/vipster.icns
RC_ICONS = resources/vipster.ico

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
        --embed-file $$PWD/resources/atom.frag@atom.frag \
        --embed-file $$PWD/resources/atom.vert@atom.vert \
        --embed-file $$PWD/resources/bond.frag@bond.frag \
        --embed-file $$PWD/resources/bond.vert@bond.vert \
        --embed-file $$PWD/resources/cell.frag@cell.frag \
        --embed-file $$PWD/resources/cell.vert@cell.vert

    copydata.files = $$PWD/vipster.html\
                     $$PWD/vipster.css\
                     $$PWD/vipster_setup.js
    copydata.commands = $(COPY_DIR) $$copydata.files $$OUT_PWD
    first.depends = $(first) copydata
    export(first.depends)
    export(copydata.commands)
    QMAKE_EXTRA_TARGETS += first copydata

    DISTFILES += \
        vipster.html \
        vipster_setup.js \
        vipster.css

} else {
    QT += core gui widgets
    TARGET = vipster

    SOURCES += main.cpp\
        mainwindow.cpp \
        glwidget.cpp \
        guiwrapper.cpp \
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
} else:win32:CONFIG(debug, debug|release) {
    LIBS += -L$$OUT_PWD/../libvipster/debug/ -lvipster
} else:unix|wasm {
    LIBS += -L$$OUT_PWD/../libvipster/ -lvipster
}

unix:!macx:CONFIG(debug, debug|release) {
    LIBS += -lgcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage -O0
    QMAKE_LFLAGS += -fprofile-arcs -ftest-coverage -O0
}
