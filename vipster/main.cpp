#include "mainwindow.h"
#include <QApplication>
#include <QCommandLineParser>
#include <iostream>
#include <iowrapper.h>

int main(int argc, char *argv[])
{
    bool parseFile = false;
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
    for(auto &kv: Vipster::IOPlugins)
    {
        const char* fmt = kv.second->argument.c_str();
        if(p.isSet(fmt)){
            parseFile = true;
            MainWindow w(Vipster::readFile(p.values(fmt).at(0).toStdString(),kv.first)->mol);
            w.show();
            return a.exec();
        }
    }
    if(!parseFile){
        MainWindow w;
        w.show();
        return a.exec();
    }
}
