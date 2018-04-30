#include "mainwindow.h"
#include "iowrapper.h"
#include "CLI11.hpp"
#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>
//TODO
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
    QApplication qapp(argc, argv);
    QApplication::setApplicationName("Vipster");
    QApplication::setApplicationVersion("1.9a");

    // main parser + data-targets
    CLI::App app{"Vipster " + QApplication::applicationVersion().toStdString()};
    // TODO
//    app.allow_extras(true);
    std::map<IOFmt, std::vector<std::string>> fmt_files{};
    std::map<CLI::Option*, IOFmt> fmt_opts{};
    for(auto& fmt: IOPlugins){
        // parser
        auto opt = app.add_option("--" + fmt.second->command,
                                  fmt_files[fmt.first],
                                  fmt.second->name);
        fmt_opts[opt] = fmt.first;
        opt->group("Parse files");
        opt->check(CLI::ExistingFile);
    }

    // conversion parser + data + options
    auto convert = app.add_subcommand("convert", "Directly convert a file");
    std::string fmt_group{"Available formats (r: parsing, w: writing)"};
    struct{
        std::vector<std::string> input;
        std::vector<std::string> output;
        std::vector<std::string> kpoints;
        std::string param;
        // TODO
//        std::unique_ptr<BaseConfig> config{};
    }conv_data;
    auto conv_kpoints = convert->add_option("-k", conv_data.kpoints,
        "Specify k-points (defaults to parsed mesh or gamma-point)");

    auto conv_param = convert->add_option("-p", conv_data.param,
        "Specify parameter set (defaults to parsed one, if present)");
    conv_param->expected(1);

    auto conv_in = convert->add_option("in", conv_data.input,
                                       "fmt and filename of input");
    conv_in->required(true);
    conv_in->expected(2);
    auto conv_out = convert->add_option("out", conv_data.output,
                                        "fmt and filename of output");
    conv_out->required(true);
    conv_out->expected(2);

    // do the parsing
    CLI11_PARSE(app, argc, argv);

    // execution
    if(app.got_subcommand(convert)){
        // TODO: as callback! should catch the throws
        IOFmt fmt_in, fmt_out;
        auto pos_in = std::find_if(IOPlugins.begin(), IOPlugins.end(),
                                   [&](const decltype(IOPlugins)::value_type& pl){
                                       return pl.second->command ==
                                               conv_data.input[0];
                                   });
        if(pos_in == IOPlugins.end()){
            throw CLI::ParseError("Invalid input format: "+conv_data.input[0], 1);
        }else{
            fmt_in = pos_in->first;
        }
        auto pos_out = std::find_if(IOPlugins.begin(), IOPlugins.end(),
                                    [&](const decltype(IOPlugins)::value_type& pl){
                                        return (pl.second->command ==
                                                conv_data.output[0]) &&
                                               (pl.second->writer != nullptr);
                                    });
        if(pos_out == IOPlugins.end()){
            throw CLI::ParseError("Invalid output format: "+conv_data.output[0], 1);
        }else{
            fmt_out = pos_out->first;
        }
        auto [mol, fmt, param] = readFile(conv_data.input[1], fmt_in);
        if(!conv_data.param.empty()){
            auto tmp = params.equal_range(fmt_in);
            auto pos_par = std::find_if(tmp.first, tmp.second,
                                        [&](const decltype(params)::value_type& p){
                                            return p.second->name == conv_data.param;
                                        });
            if(pos_par == params.end()){
                throw CLI::ParseError("Invalid parameter \""+conv_data.param+
                                      "\" for format "+conv_data.input[0], 1);
            }
            param = pos_par->second->copy();
        }
        if(!conv_data.kpoints.empty()){
            const auto& kpoints = conv_data.kpoints;
            const std::string kp_err = "KPoints should be one of:\n"
                                       "gamma\n"
                                       "mpg x y z sx sy sz\n"
                                       "";
            if(kpoints[0] == "gamma"){
                if(kpoints.size() > 1){
                    throw CLI::ParseError("gamma takes no arguments\n\n"+kp_err, 1);
                }
                mol.setKPoints(KPoints{KPointFmt::Gamma, {}, {}});
            }else if(kpoints[0] == "mpg"){
                const std::string mpg_err = "mpg takes three integer and three "
                                            "float arguments\n\n";
                if(kpoints.size() != 7){
                    throw CLI::ParseError(mpg_err+kp_err, 1);
                }
                try {
                    mol.setKPoints(KPoints{KPointFmt::MPG, {
                                       std::stoi(kpoints[1]),
                                       std::stoi(kpoints[2]),
                                       std::stoi(kpoints[3]),
                                       std::stof(kpoints[4]),
                                       std::stof(kpoints[5]),
                                       std::stof(kpoints[6]),
                                   }, {}});
                } catch (...) {
                    throw CLI::ParseError(mpg_err+kp_err, 1);
                }{}
            //TODO: discrete
            }else{
                throw CLI::ParseError("Invalid KPoint style\n"+kp_err, 1);
            }
        }
        writeFile(conv_data.output[1], fmt_out, mol, param.get());
    }else{
        std::vector<IO::Data> data{};
        // TODO: try to guess formats from app.remaining()
        // parse files
        for(auto& op_fmt: fmt_opts){
            for(auto& fn: op_fmt.first->results()){
                data.push_back(readFile(fn, op_fmt.second));
            }
        }
        // launch GUI
        if(data.size()){
            MainWindow w{std::move(data)};
            w.show();
            return QApplication::exec();
        }else{
            MainWindow w{};
            w.show();
            return QApplication::exec();
        }
    }
}
