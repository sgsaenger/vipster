#include "version.h"
#include "mainwindow.h"
#include "vipster/fileio.h"
#include "vipster/configfile.h"
#include "CLI/CLI.hpp"
#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>
#include <QOpenGLContext>
#include <fmt/format.h>

#ifdef USE_PYTHON
#pragma push_macro("slots")
#undef slots
#include <pybind11/embed.h>
#include "vipster/global.py.h"
#pragma pop_macro("slots")
#endif

#ifdef USE_LAMMPS
#include "toolwidgets/lammpswidget_aux/run.lmp.h"
#endif

using namespace Vipster;

// setup and launch GUI
[[noreturn]] void launchVipster(int argc, char *argv[],
                                std::vector<IOTuple>&& data,
                                ConfigState&& state){
    QApplication qapp(argc, argv);
    QApplication::setApplicationName("Vipster");
    QApplication::setApplicationVersion(VIPSTER_VERSION);
    std::cout << "Vipster v" VIPSTER_VERSION << std::endl;
    QSurfaceFormat format = QSurfaceFormat::defaultFormat();
    if(QOpenGLContext::openGLModuleType() == QOpenGLContext::LibGL){
        // try to get a 3.3core context on desktop
        format.setVersion(3,3);
        format.setProfile(QSurfaceFormat::CoreProfile);
        format.setRenderableType(QSurfaceFormat::OpenGL);
    }else{
        // or an es3.0 context on mobile
        format.setVersion(3,0);
        format.setRenderableType(QSurfaceFormat::OpenGLES);
    }
    format.setAlphaBufferSize(8);
    format.setSamples(8);
    QSurfaceFormat::setDefaultFormat(format);
    QObject::connect(&qapp, &QApplication::aboutToQuit, &qapp, [&](){saveConfig(state);});

    MainWindow w{QDir::currentPath(), state, std::move(data)};
    w.show();
    throw CLI::RuntimeError(QApplication::exec());
}

// register `convert` subcommand
void addSubcommandConvert(CLI::App& app, const ConfigState& state){
    // create a CLI11 subcommand
    auto convert = app.add_subcommand("convert", "Convert a file from one format to another");
    convert->set_help_all_flag();

    // storage for parsed options
    static struct{
        std::string in_fmt;
        std::string in_fn;
        std::string out_fmt;
        std::string out_fn;
        std::vector<std::string> kpoints;
        std::string param;
        std::string preset;
    }conv_data;

    // aliases for config state members
    const PluginList    &plugins = std::get<2>(state);
    const ParameterMap  &params = std::get<3>(state);
    const PresetMap     &presets = std::get<4>(state);

    // formats
    convert->add_flag("--list-fmt",
                      [&](size_t){
                          std::cout << "Available formats (r: parsing, w: writing)\n\n";
                          for(const auto& plug: plugins){
                              std::cout << plug->command << "\t("
                                        << (plug->parser?'r':' ')
                                        << (plug->writer?'w':' ')
                                        << ")\t" << plug->name << '\n';
                          }
                          throw CLI::Success();
                      },
                      "Display available formats");
    // K-points
    constexpr const char* kp_err = "KPoints should be one of:\n\n"
        "gamma\t\t\t\tGamma-point only\n"
        "mpg x y z sx sy sz\t\tMonkhorst-pack grid of size x*y*z with offset (sx,sy,sz)\n"
        "disc N B C [x y z w](Nx)\tDiscrete grid with N points, each given with position (x,y,z) and weight w.\n"
        "\t\t\t\tB: Toggle Band-mode (0,1); C: Toggle crystal-coordinates (0,1)"
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
                          std::cout << Vipster::ParametersAbout << std::endl;
                          throw CLI::Success();
                      },
                      "Display help for parameter sets");
    convert->add_flag("--list-param",
                      [&](size_t){
                          auto printFmt = [&](const Plugin* fmt){
                              for(const auto& pair: params.at(fmt)){
                                  std::cout << pair.first << '\n';
                              }
                          };
                          for(const auto& plug: plugins){
                              if(plug->makeParam){
                                  std::cout << plug->command << ": "
                                            << plug->name << "\n";
                                  printFmt(plug);
                                  std::cout << '\n';
                              }
                          }
                          throw CLI::Success();
                      },
                      "List available parameter sets");

    // IO-presets
    convert->add_option("-c,--preset", conv_data.preset,
                        "Specify behavior-preset for output plugin");
    convert->add_flag("--help-preset",
                      [](size_t){
                          std::cout << Vipster::PresetsAbout << std::endl;
                          throw CLI::Success();
                      },
                      "Display help for output-behavior-presets");
    convert->add_flag("--list-preset",
                      [&](size_t){
                          auto printFmt = [&](const Plugin* fmt){
                              for(const auto& pair: presets.at(fmt)){
                                  std::cout << pair.first << '\n';
                              }
                          };
                          for(const auto& plug: plugins){
                              if(plug->makePreset){
                                  std::cout << plug->command << ": "
                                            << plug->name << "\n";
                                  printFmt(plug);
                                  std::cout << '\n';
                              }
                          }
                          throw CLI::Success();
                      },
                      "List available output-behavior-presets");

    // main file arguments
    convert->add_option("in_fmt", conv_data.in_fmt,
                        "format of input file")->required(true);
    convert->add_option("in_fn", conv_data.in_fn,
                        "input filename")->required(true);
    convert->add_option("out_fmt", conv_data.out_fmt,
                        "format of output file")->required(true);
    convert->add_option("out_fn", conv_data.out_fn,
                        "output filename")->required(true);

    convert->callback([&](){
        // determine/check in&out formats
        const Plugin *fmt_in, *fmt_out;
        const auto IOCmdIn = [&](){
            std::map<std::string, const Plugin*> fmts_in;
            for(const auto& plug: plugins){
                if(plug->parser){
                    fmts_in[plug->command] = plug;
                }
            }
            return fmts_in;
        }();
        const auto IOCmdOut = [&](){
            std::map<std::string, const Plugin*> fmts_out;
            for(const auto& plug: plugins){
                if(plug->writer){
                    fmts_out[plug->command] = plug;
                }
            }
            return fmts_out;
        }();
        auto pos_in = IOCmdIn.find(conv_data.in_fmt);
        if(pos_in == IOCmdIn.end()){
            throw CLI::ParseError("Invalid input format: "+conv_data.in_fmt, 1);
        }else{
            fmt_in = pos_in->second;
        }
        auto pos_out = IOCmdOut.find(conv_data.out_fmt);
        if(pos_out == IOCmdOut.end()){
            throw CLI::ParseError("Invalid output format: "+conv_data.out_fmt, 1);
        }else{
            fmt_out = pos_out->second;
        }
        // read input file
        auto [mol, param, data] = readFile(conv_data.in_fn, fmt_in);
        std::optional<Preset> preset{};
        // create Parameter set if requested/required
        if(fmt_out->makeParam){
            std::string par_name;
            if(!conv_data.param.empty()){
                par_name = conv_data.param;
            }else if(!param){
                par_name = "default";
            }
            if(!par_name.empty()){
                const auto& pos = params.at(fmt_out).find(par_name);
                if(pos == params.at(fmt_out).end()){
                    throw CLI::ParseError("Invalid parameter \""+par_name+
                                          "\" for format "+conv_data.out_fmt, 1);
                }
                param = pos->second;
            }
        }
        // create IO preset set if requested/required
        if(fmt_out->makePreset){
            std::string pres_name;
            if(!conv_data.preset.empty()){
                pres_name = conv_data.preset;
            }else{
                pres_name = "default";
            }
            auto pos = presets.at(fmt_out).find(pres_name);
            if(pos == presets.at(fmt_out).end()){
                throw CLI::ParseError("Invalid IO preset \""+pres_name+
                                      "\" for format "+conv_data.out_fmt, 1);
            }
            preset = pos->second;
        }
        // parse k-points
        if(!conv_data.kpoints.empty()){
            const auto& kpoints = conv_data.kpoints;
            if(kpoints[0] == "gamma"){
                mol.kpoints.active = KPoints::Fmt::Gamma;
                if(kpoints.size() > 1){
                    throw CLI::ParseError("K-Point \"gamma\" takes no arguments", 1);
                }
            }else if(kpoints[0] == "mpg"){
                mol.kpoints.active = KPoints::Fmt::MPG;
                const std::string mpg_err = "Monkhorst-Pack grid expects three integer "
                                            "and three float arguments";
                if(kpoints.size() != 7){
                    throw CLI::ParseError(mpg_err, 1);
                }
                try {
                    mol.kpoints.mpg = {std::stoi(kpoints[1]),
                                       std::stoi(kpoints[2]),
                                       std::stoi(kpoints[3]),
                                       std::stod(kpoints[4]),
                                       std::stod(kpoints[5]),
                                       std::stod(kpoints[6])};
                } catch (...) {
                    throw CLI::ParseError(mpg_err, 1);
                }
            }else if(kpoints[0] == "disc"){
                mol.kpoints.active = KPoints::Fmt::Discrete;
                if(kpoints.size() < 4){
                    throw CLI::ParseError("K-Point \"disc\" takes at least three arguments", 1);
                }
                auto N = static_cast<size_t>(std::stoi(kpoints[1]));
                if(kpoints.size() != 4+4*N){
                    throw CLI::ParseError("K-Point \"disc\" takes exactly three + 4*N arguments", 1);
                }
                auto &properties = mol.kpoints.discrete.properties;
                if(std::stoi(kpoints[2])){
                    properties |= KPoints::Discrete::band;
                }
                if(std::stoi(kpoints[3])){
                    properties |= KPoints::Discrete::crystal;
                }
                auto& list = mol.kpoints.discrete.kpoints;
                list.resize(N);
                size_t i = 4;
                for(auto& kpoint: list){
                    kpoint = {{std::stod(kpoints[i]),
                               std::stod(kpoints[i+1]),
                               std::stod(kpoints[i+2])},
                              std::stod(kpoints[i+3]),
                             };
                    i+=4;
                }
            }else{
                throw CLI::ParseError(std::string{"Invalid KPoint style\n"}+kp_err, 1);
            }
        }
        writeFile(conv_data.out_fn, fmt_out, mol, std::nullopt, param, preset);
        throw CLI::Success();
    });

}

#ifdef USE_PYTHON
// create embedded python module
PYBIND11_EMBEDDED_MODULE(vipster, m) {}
#endif

// command-line handling
int main(int argc, char *argv[])
{
    // read user-defined settings and make state known
    auto state = Vipster::readConfig();
    const PluginList    &plugins = std::get<2>(state);
    const ParameterMap  &params = std::get<3>(state);
    const PresetMap     &presets = std::get<4>(state);

#ifdef USE_PYTHON
    // instance the python-interpreter, keep it alive for the program's duration
    pybind11::scoped_interpreter interp{};
    auto vipster_module = py::module::import("vipster");
    // register types and state in module
    Py::setupVipster(vipster_module, state, false);
#endif

    // main parser + data-targets
    CLI::App app{"Vipster v" VIPSTER_VERSION};
    app.set_help_all_flag("-H,--help-all", "Print help for this and all subcommands");
    app.allow_extras(true);
    std::map<const Plugin*, std::vector<std::string>> fmt_files{};
    std::map<CLI::Option*, const Plugin*> fmt_opts{};
    for(auto& fmt: plugins){
        // parser
        if(!fmt->parser) continue;
        try{
            auto opt = app.add_option("--" + fmt->command,
                                      fmt_files[fmt],
                                      fmt->name);
            fmt_opts[opt] = fmt;
            opt->group("Parse file(s) or stdin ('-')");
        }catch(CLI::OptionAlreadyAdded &e){
            std::cerr << fmt::format("Unable to activate plugin {}: {}",
                                     fmt->name, e.what()) << std::endl;
        }
    }
    app.callback([&](){
        if(!app.get_subcommands().empty()){
            return;
        }
        std::vector<IOTuple> data{};
        if(app.remaining_size()!=0){
            for(const auto& file: app.remaining()){
                if(file == "-"){
                    std::cout << "Cannot parse stdin without explicit format" << std::endl;
                    throw CLI::RuntimeError{1};
                }
                if(const auto plug = guessFmt(file, std::get<2>(state))){
                    data.push_back(readFile(file, plug));
                }else{
                    std::cout << fmt::format("Could not deduce format of file \"{}\""
                                             "\nPlease specify format explicitely", file)
                              << std::endl;
                    throw CLI::RuntimeError{1};
                }
            }
        }
        // parse files
        for(auto& op_fmt: fmt_opts){
            for(const auto& fn: op_fmt.first->results()){
                try{
                    data.push_back(readFile(fn, op_fmt.second));
                }catch(const Vipster::IOError &e){
                    std::cout << e.what() << std::endl;
                    throw CLI::RuntimeError{1};
                }
            }
        }
        // launch GUI
        launchVipster(argc, argv, std::move(data), std::move(state));
    });

    // add `convert` subcommand
    addSubcommandConvert(app, state);

#if defined(USE_LAMMPS) && defined(USE_MPI)
    auto lmp = app.add_subcommand("lammps_mpi_slave");
    lmp->group("");
    lmp->callback([&](){
        // launch lammps mpi slaves
        Lammps::runSlave();
        throw CLI::Success();
    });
#endif

    // do the parsing
    CLI11_PARSE(app, argc, argv)
}
