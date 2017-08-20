#include "mainwindow.h"
#include <QApplication>
#include <QCommandLineParser>
#include <iostream>
#include <iowrapper.h>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    a.setApplicationName("Vipster");
    a.setApplicationVersion("1.9a");
    QCommandLineParser p;
    p.setApplicationDescription("Vipster");
    p.addHelpOption();
    p.addVersionOption();
    for(auto &kv: Vipster::IOPlugins)
    {
        p.addOption({QString::fromStdString(kv.second->argument),
                     QString::fromStdString(kv.second->name),
                     "files"});
    }
    p.process(a);
    if(p.isSet("xyz")){
        MainWindow w(Vipster::readFile(p.values("xyz").at(0).toStdString(),Vipster::IOFmt::XYZ)->mol);
        w.show();
        return a.exec();
    }else if(p.isSet("pwi")){
        MainWindow w(Vipster::readFile(p.values("pwi").at(0).toStdString(),Vipster::IOFmt::PWI)->mol);
        w.show();
        return a.exec();
    }else{
        MainWindow w;
        w.show();
        return a.exec();
    }
}
