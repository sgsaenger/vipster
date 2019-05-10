#include "version.h"
#include "mainwindow.h"
#include "io.h"
#include "configfile.h"
#include "CLI11.hpp"
#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>

using namespace Vipster;

[[noreturn]] void launchVipster(int argc, char *argv[], std::vector<IO::Data>&& data){
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setSamples(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(format);
    QApplication qapp(argc, argv);
    QApplication::setApplicationName("Vipster");
    QObject::connect(&qapp, &QApplication::aboutToQuit, &qapp, [](){saveConfig();});
    if(!data.empty()){
        MainWindow w{QDir::currentPath(), std::move(data)};
        w.show();
        throw CLI::RuntimeError(QApplication::exec());
    }else{
        MainWindow w{QDir::currentPath()};
        w.show();
        throw CLI::RuntimeError(QApplication::exec());
    }
}

int main(int argc, char *argv[])
{
    // main parser + data-targets
    Vipster::readConfig();
    CLI::App app{"Vipster v" VIPSTER_VERSION "b"};
    app.allow_extras(true);
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
    app.callback([&](){
        if(!app.get_subcommands().empty()){
            return;
        }
        std::vector<IO::Data> data{};
        if(app.remaining_size()!=0){
            const auto IOExt = [](){
                std::map<std::string, IOFmt> ext;
                for(const auto& pair: IOPlugins){
                    ext[pair.second->extension] = pair.first;
                }
                return ext;
            }();
            for(const auto& file: app.remaining()){
                auto pos = file.find_last_of('.');
                if(pos != file.npos){
                    std::string ext = file.substr(pos+1);
                    auto pos = IOExt.find(ext);
                    if(pos != IOExt.end()){
                        data.push_back(readFile(file, pos->second));
                    }else{
                        throw CLI::ParseError(
                                    "Could not deduce format of file "+file+
                                    "\nPlease specify format explicitely", 1);
                    }
                }else{
                    throw CLI::ParseError(
                                "Could not deduce format of file "+file+
                                "\nPlease specify format explicitely", 1);
                }
            }
        }
        // parse files
        for(auto& op_fmt: fmt_opts){
            for(const auto& fn: op_fmt.first->results()){
                data.push_back(readFile(fn, op_fmt.second));
            }
        }
        // launch GUI
        launchVipster(argc, argv, std::move(data));
    });

    // conversion parser + data + options
    auto convert = app.add_subcommand("convert", "Directly convert a file");
    struct{
        std::vector<std::string> input;
        std::vector<std::string> output;
        std::vector<std::string> kpoints;
        std::string param;
        std::string config;
    }conv_data;

    // formats
    convert->add_flag("--list-fmt",
                      [](size_t){
                          std::cout << "Available formats (r: parsing, w: writing)\n\n";
                          for(const auto& pair: IOPlugins){
                              std::cout << pair.second->command << "\t(r"
                                        << ((pair.second->writer!=nullptr)?'w':' ')
                                        << ")\t" << pair.second->name << '\n';
                          }
                          throw CLI::Success();
                      },
                      "Display available formats");
    // K-points
    constexpr const char* kp_err = "KPoints should be one of:\n\n"
        "gamma\t\t\tGamma-point only\n"
        "mpg x y z sx sy sz\t\tMonkhorst-pack grid of size x*y*z with offset (sx,sy,sz)\n"
        "disc N B C [x y z w](Nx)\tDiscrete grid with N points, each given with position (x,y,z) and weight w.\n"
        "\t\t\tB: Toggle Band-mode (0,1); C: Toggle crystal-coordinates (0,1)"
        "";
    convert->add_option("-k,--kpoints", conv_data.kpoints,
        "Specify k-points (defaults to parsed mesh or gamma-point)");
    convert->add_flag("--help-kpoints",
                      [](size_t){
                          std::cout << kp_err << std::endl;
                          throw CLI::Success();
                      },
                      "Display help for k-point specification");

    // parameter sets
    convert->add_option("-p,--param", conv_data.param,
                        "Specify parameter set (defaults to parsed one, if present)");
    convert->add_flag("--help-param",
                      [](size_t){
                          std::cout << Vipster::IO::ParametersAbout << std::endl;
                          throw CLI::Success();
                      },
                      "Display help for parameter sets");
    convert->add_flag("--list-param",
                      [](size_t){
                          auto printFmt = [](IOFmt fmt){
                              for(const auto& pair: params[fmt]){
                                  std::cout << pair.first << '\n';
                              }
                          };
                          for(const auto& pair: IOPlugins){
                              if(pair.second->arguments&IO::Plugin::Param){
                                  std::cout << pair.second->command << ": "
                                            << pair.second->name << "\n";
                                  printFmt(pair.first);
                                  std::cout << '\n';
                              }
                          }
                          throw CLI::Success();
                      },
                      "List available parameter sets");

    // IO-configs
    convert->add_option("-c,--config", conv_data.config,
                        "Specify behavior-preset for output plugin");
    convert->add_flag("--help-config",
                      [](size_t){
                          std::cout << Vipster::IO::ConfigsAbout << std::endl;
                          throw CLI::Success();
                      },
                      "Display help for output-behavior-presets");
    convert->add_flag("--list-config",
                      [](size_t){
                          auto printFmt = [](IOFmt fmt){
                              for(const auto& pair: configs[fmt]){
                                  std::cout << pair.first << '\n';
                              }
                          };
                          for(const auto& pair: IOPlugins){
                              if(pair.second->arguments&IO::Plugin::Config){
                                  std::cout << pair.second->command << ": "
                                            << pair.second->name << "\n";
                                  printFmt(pair.first);
                                  std::cout << '\n';
                              }
                          }
                          throw CLI::Success();
                      },
                      "List available output-behavior-presets");

    // main arguments
    auto conv_in = convert->add_option("in", conv_data.input,
                                       "fmt and filename of input");
    conv_in->required(true);
    conv_in->expected(2);
    auto conv_out = convert->add_option("out", conv_data.output,
                                        "fmt and filename of output");
    conv_out->required(true);
    conv_out->expected(2);

    convert->callback([&](){
        // determine/check in&out formats
        IOFmt fmt_in, fmt_out;
        const auto IOCmdIn = [](){
            std::map<std::string, IOFmt> fmts_in;
            for(const auto& pair: IOPlugins){
                fmts_in[pair.second->command] = pair.first;
            }
            return fmts_in;
        }();
        const auto IOCmdOut = [](){
            std::map<std::string, IOFmt> fmts_out;
            for(const auto& pair: IOPlugins){
                if(pair.second->writer!=nullptr){
                    fmts_out[pair.second->command] = pair.first;
                }
            }
            return fmts_out;
        }();
        auto pos_in = IOCmdIn.find(conv_data.input[0]);
        if(pos_in == IOCmdIn.end()){
            throw CLI::ParseError("Invalid input format: "+conv_data.input[0], 1);
        }else{
            fmt_in = pos_in->second;
        }
        auto pos_out = IOCmdOut.find(conv_data.output[0]);
        if(pos_out == IOCmdOut.end()){
            throw CLI::ParseError("Invalid output format: "+conv_data.output[0], 1);
        }else{
            fmt_out = pos_out->second;
        }
        // read input
        auto [mol, param, data] = readFile(conv_data.input[1], fmt_in);
        std::unique_ptr<IO::BaseConfig> config{nullptr};
        auto arguments = IOPlugins.at(fmt_out)->arguments;
        if(arguments & IO::Plugin::Args::Param){
            std::string par_name;
            if(!conv_data.param.empty()){
                par_name = conv_data.param;
            }else if(!param){
                par_name = "default";
            }
            if(!par_name.empty()){
                const auto& pos = params[fmt_out].find(par_name);
                if(pos == params[fmt_out].end()){
                    throw CLI::ParseError("Invalid parameter \""+par_name+
                                          "\" for format "+conv_data.output[0], 1);
                }
                param = pos->second->copy();
            }
        }
        if(arguments & IO::Plugin::Args::Config){
            std::string conf_name;
            if(!conv_data.config.empty()){
                conf_name = conv_data.config;
            }else{
                conf_name = "default";
            }
            auto pos = configs[fmt_out].find(conf_name);
            if(pos == configs[fmt_out].end()){
                throw CLI::ParseError("Invalid configuration preset \""+conf_name+
                                      "\" for format "+conv_data.output[0], 1);
            }
            config= pos->second->copy();
        }
        if(!conv_data.kpoints.empty()){
            const auto& kpoints = conv_data.kpoints;
            if(kpoints[0] == "gamma"){
                if(kpoints.size() > 1){
                    throw CLI::ParseError("K-Point \"gamma\" takes no arguments", 1);
                }
                mol.setKPoints(KPoints{KPoints::Fmt::Gamma, {}, {}});
            }else if(kpoints[0] == "mpg"){
                const std::string mpg_err = "Monkhorst-Pack grid expects three integer "
                                            "and three float arguments";
                if(kpoints.size() != 7){
                    throw CLI::ParseError(mpg_err, 1);
                }
                try {
                    mol.setKPoints(KPoints{KPoints::Fmt::MPG, {
                                       std::stoi(kpoints[1]),
                                       std::stoi(kpoints[2]),
                                       std::stoi(kpoints[3]),
                                       std::stof(kpoints[4]),
                                       std::stof(kpoints[5]),
                                       std::stof(kpoints[6]),
                                   }, {}});
                } catch (...) {
                    throw CLI::ParseError(mpg_err, 1);
                }
            //TODO: discrete
            }else if(kpoints[0] == "disc"){
                if(kpoints.size() < 4){
                    throw CLI::ParseError("K-Point \"disc\" takes at least three arguments", 1);
                }
                auto N = static_cast<size_t>(std::stoi(kpoints[1]));
                if(kpoints.size() != 4+4*N){
                    throw CLI::ParseError("K-Point \"disc\" takes exactly three + 4*N arguments", 1);
                }
                uint8_t properties{};
                if(std::stoi(kpoints[2])){
                    properties |= KPoints::Discrete::band;
                }
                if(std::stoi(kpoints[3])){
                    properties |= KPoints::Discrete::crystal;
                }
                mol.setKPoints({KPoints::Fmt::Discrete, {}, {properties}});
                auto& list = mol.getKPoints().discrete.kpoints;
                list.resize(N);
                size_t i = 4;
                for(auto& kpoint: list){
                    kpoint = {{std::stof(kpoints[i]),
                               std::stof(kpoints[i+1]),
                               std::stof(kpoints[i+2])},
                              std::stof(kpoints[i+3]),
                             };
                    i+=4;
                }
            }else{
                throw CLI::ParseError(std::string{"Invalid KPoint style\n"}+kp_err, 1);
            }
        }
        writeFile(conv_data.output[1], fmt_out, mol, param.get(), config.get());
        throw CLI::Success();
    });

    // do the parsing
    CLI11_PARSE(app, argc, argv);
}
