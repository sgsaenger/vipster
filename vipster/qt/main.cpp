#include "mainwindow.h"
#include "iowrapper.h"
#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>
#include <iostream>

using namespace Vipster;

int main(int argc, char *argv[])
{
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setSamples(8);
    format.setAlphaBufferSize(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(format);
    QApplication a(argc, argv);
    a.setApplicationName("Vipster");
    a.setApplicationVersion("1.8a");
    QCommandLineParser p;
    p.setApplicationDescription(
        "Vipster, a graphical molecular structure viewer and editor\n\n"
        "Load files via:\n"
        "--<fmt> [files...]\n\n"
//        "Convert files via:\n"
//        "--convert --<fmt_in> <file_in> --<fmt_out> <file_out>"
    );
    p.addHelpOption();
    p.addVersionOption();
    // Add conversion options
    // TODO:
//    p.addOption({{"convert","c"}, "Convert from one format to another"});
//    p.addOption({{"showWriters","w"}, "List available writers for conversion"});
    /*
     * handle optional arguments for conversion: Kpoints, parameters, (config?)
     */
    // Add parser options
    QList<QCommandLineOption> fmt_options{};
    for(auto &kv: IOPlugins)
    {
        fmt_options.push_back({QString::fromStdString(kv.second->command),
                               QString::fromStdString(kv.second->name)});
    }
    p.addOptions(fmt_options);
    // Process arguments:
    p.process(a);
    if(!p.optionNames().size()){
        // No options -> launch empty GUI
        if(p.positionalArguments().size()){
            p.showHelp(1);
        }
        MainWindow w{};
        w.show();
        return a.exec();
    } else {
        if(p.isSet("w")){
            std::cout << "Vipster\n\nAvailable parsers:\n";
            for(auto& p:IOPlugins){
                if(p.second->writer){
                    std::cout << "--" << p.second->command << '\t' << p.second->name << '\n';
                }
            }
            return 0;
        }
        else if(p.isSet("convert")){
            //TODO: perform conversion!
            return 1;
        } else {
            if(p.optionNames().size() != 1){
                p.showHelp(1);
            }
            std::vector<IO::Data> data;
            IOFmt fmt;
            for(auto& kv:IOPlugins){
                if(p.isSet(kv.second->command.c_str())){
                    fmt = kv.first;
                }
            }
            for(auto& file:p.positionalArguments()){
                data.push_back(readFile(file.toStdString(), fmt));
            }
            MainWindow w{std::move(data)};
            w.show();
            return a.exec();
        }
    }
}
