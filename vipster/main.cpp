#include "mainwindow.h"
#include "iowrapper.h"
#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>
#include <iostream>

int main(int argc, char *argv[])
{
    bool parseFile = false;
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setSamples(8);
    format.setAlphaBufferSize(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(format);
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
