#include "mainwindow.h"
#include <QApplication>
#include <QCommandLineParser>
#include "libvipster.h"
#include <iostream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    a.setApplicationName("Vipster");
    a.setApplicationVersion("0.9");
    QCommandLineParser p;
    p.setApplicationDescription("Vipster");
    p.addHelpOption();
    p.addVersionOption();
    p.addOption({"xyz","Parse xyz <files>.","files"});
    p.process(a);
    Vipster::Molecule m;
    if(p.isSet("xyz")){
        std::tie(m,std::ignore) = Vipster::readFile(p.values("xyz").at(0).toStdString(),"xyz");
        MainWindow w(m);
        w.show();
        return a.exec();
    }else{
        MainWindow w;
        w.show();
        return a.exec();
    }
}
