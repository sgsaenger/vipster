#include "lammpswidget.lmp.h"
#include "ui_lammpswidget.lmp.h"
#include "../mainwindow.h"
#include "lammpswidget_aux/run.lmp.h"

#include "vipster/plugins/lmpinput.h"

#include <thread>

#include <fmt/format.h>

#include <QMessageBox>

#include "lammps.h"
#include "info.h"
#include "exceptions.h"
#include "vipsterapplication.h"

using namespace Vipster;
using namespace Vipster::Lammps;
using namespace LAMMPS_NS;
namespace fs = std::filesystem;

LammpsWidget::LammpsWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::LammpsWidget)
{
    ui->setupUi(this);
    ui->customFrame->setHidden(true);
    // set validators for minimize settings
    ui->etolInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->ftolInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->maxItInput->setValidator(new QIntValidator{1, 10000000});
    ui->maxevalInput->setValidator(new QIntValidator{1, 10000000});
    // set validators for MD settings
    ui->tempInput->setValidator(new QDoubleValidator{0.1, 100000, 5});
    ui->pressInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->stepInput->setValidator(new QIntValidator{1, 10000000});
    ui->timeInput->setValidator(new QDoubleValidator{0.00001, 1000, 5});
    // enable parallelization as available
    auto nproc = std::max(1u, std::thread::hardware_concurrency());
//#ifdef USE_MPI
    ui->MPISpin->setDisabled(true);
//#else
//    ui->MPISpin->setEnabled(true);
//    ui->MPISpin->setMinimum(1);
//    ui->MPISpin->setMaximum(nproc);
//#endif
    ui->OMPSpin->setDisabled(true);
    ui->OMPSpin->setMaximum(nproc);
    ui->GPUSpin->setDisabled(true);
    int i=0;
    while(LAMMPS::installed_packages[i] != nullptr){
        if(std::string{"USER-OMP"} == LAMMPS::installed_packages[i]){
            ui->OMPSpin->setEnabled(true);
        }else if(std::string{"GPU"} == LAMMPS::installed_packages[i]){
            ui->GPUSpin->setEnabled(true);
        }
        ++i;
    }

    int flag;
    MPI_Initialized(&flag);
    if(!flag){
        int argc{0};
        char **argv{nullptr};
#ifdef USE_MPI
        int level{0};
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
        if(level < MPI_THREAD_MULTIPLE){
            ui->MPISpin->setDisabled(true);
        }
#else
        // LAMMPS-MPI stubs doesn't implement the threaded functions
        MPI_Init(&argc, &argv);
#endif
    }

    // register installed FF styles
    std::array<char*, 5> lmparg{
        nullptr,
        "-screen", "none",
        "-log", "none"
    };
    LAMMPS lmp{5, lmparg.data(), MPI_COMM_WORLD};
    Info info{&lmp};
    for(const auto& style: info.get_available_styles("pair")){
        ui->pairSel->addItem(style.c_str());
    }
    for(const auto& style: info.get_available_styles("bond")){
        ui->bondSel->addItem(style.c_str());
    }
    for(const auto& style: info.get_available_styles("angle")){
        ui->angleSel->addItem(style.c_str());
    }
    for(const auto& style: info.get_available_styles("dihedral")){
        ui->dihedSel->addItem(style.c_str());
    }
    for(const auto& style: info.get_available_styles("improper")){
        ui->impropSel->addItem(style.c_str());
    }
    for(const auto& style: info.get_available_styles("kspace")){
        ui->kspaceSel->addItem(style.c_str());
    }

    // register min/md styles
    for(const auto& style: info.get_available_styles("minimize")){
        ui->minSel->addItem(style.c_str());
    }
    if(info.has_style("fix", "nve"))
        ui->mdSel->addItem("nve");
    if(info.has_style("fix", "nvt"))
        ui->mdSel->addItem("nvt");
    if(info.has_style("fix", "npt"))
        ui->mdSel->addItem("npt");
    if(info.has_style("fix", "nph"))
        ui->mdSel->addItem("nph");

    // register forcefields if lammps is compatible
    for(const auto& FF: defaultForcefields()){
        // skip inserting if any style is not found in lammps
        const auto& impl = *FF.second;
        if(impl.pair.has_value()){
            // special treatment because of cutoffs
            std::string name;
            std::stringstream line{impl.pair.value()};
            bool val{true};
            while(!(line >> name).eof()){
                if(!isalpha(name[0])){
                    continue;
                }
                val &= info.has_style("pair", name);
            }
            if(!val){
                continue;
            }
        }
        auto checkStyle = [&](const std::string &type, const std::string &line){
            std::stringstream ls{line};
            std::string name;
            bool val{true};
            while(!(ls >> name).eof()){ // loop to support hybrid
                val &= info.has_style(type, name);
            }
            return val;
        };
        if(impl.bond.has_value() && !checkStyle("bond", impl.bond.value())){
            continue;
        }
        if(impl.angle.has_value() && !checkStyle("angle", impl.angle.value())){
            continue;
        }
        if(impl.dihedral.has_value() && !checkStyle("dihedral", impl.dihedral.value())){
            continue;
        }
        if(impl.improper.has_value() && !checkStyle("improper", impl.improper.value())){
            continue;
        }
        forcefields.insert(FF);
        ui->ffSel->addItem(FF.first.c_str());
    }
    ui->ffSel->addItem("Custom");
    // ensure widget is in valid state
    on_ffSel_currentIndexChanged(0);

    // TODO: register callback or create own fix
}

LammpsWidget::~LammpsWidget()
{
    MPI_Finalize();
    delete ui;
}

fs::path getLmpTmpDir()
{
    // TODO: add random token
    return getTempPath()/"lammps";
}

void LammpsWidget::on_runButton_clicked()
{
    if(!vApp.curStep().getNat()){
        return;
    }
    auto curStep = vApp.curStep().asFmt(AtomFmt::Angstrom);
    // get forcefield
    const auto FFname = ui->ffSel->currentText().toStdString();
    if(FFname == "Custom"){
        return;
    }
    const auto& FF = forcefields.at(FFname);
    // manually unset and save user locale
    std::string userLocale = setlocale(0, nullptr);
    setlocale(LC_ALL, "C");
    try{
        // prepare working directory
        auto tempdir = getLmpTmpDir();
        fs::create_directory(tempdir);
        // prepare input files
        mkGeom(curStep, *FF, tempdir);
        mkScript(curStep, *FF, tempdir);
        bool doMin = ui->calcStack->currentIndex() == 0;
        auto params = doMin ? runParams{Lammps::runParams::Mode::Min,
                                        ui->MPISpin->value(),
                                        ui->OMPSpin->value(),
                                        ui->GPUSpin->value(),
                                        ui->maxItInput->text().toULong(),
                                        ui->maxevalInput->text().toULong(),
                                        ui->etolInput->text().toDouble(),
                                        ui->ftolInput->text().toDouble()}
                            : runParams{Lammps::runParams::Mode::MD,
                                        ui->MPISpin->value(),
                                        ui->OMPSpin->value(),
                                        ui->GPUSpin->value(),
                                        ui->stepInput->text().toULong()};
        auto name = doMin ? fmt::format("(Min: {})", ui->minSel->currentText().toStdString())
                          : fmt::format("MD: {})", ui->mdSel->currentText().toStdString());
        Molecule mol{vApp.curStep(), vApp.curMol->name + name};
        auto result = runMaster(tempdir.string(), params, &mol);
        if(result.first < 0){
            QMessageBox::critical(this, "Error in LAMMPS run", QString::fromStdString(result.second)+
                                  "\nThis error may not be recoverable from.\n"
                                  "Trying to forcefully terminate LAMMPS processes, "
                                  "please restart Vipster in case of unexpected behavior or performance loss.");
        }else if(result.first > 0){
            QMessageBox::warning(this, "Error in LAMMPS run", QString::fromStdString(result.second));
        }
        // notify vipster of new mol
        vApp.newMol(std::move(mol));
    }catch(std::exception &e){
        QMessageBox::warning(this, "Error in LAMMPS run", QString::fromStdString(e.what()));
    }
    // restore user locale
    setlocale(LC_ALL, userLocale.c_str());
}

void LammpsWidget::on_helpButton_clicked()
{
    // TODO: explain this widget
}

void LammpsWidget::on_ffPrepare_clicked()
{
    const auto FFname = ui->ffSel->currentText().toStdString();
    if(FFname == "Custom"){
        return;
    }
    const auto& FF = forcefields.at(FFname);
    if(FF->prepareStep){
        try{
            auto mol = FF->prepareStep(vApp.curStep(), vApp.curMol->name);
            vApp.newMol(std::move(mol));
        }catch(const Vipster::Error &e){
            QMessageBox::warning(this, "Could not prepare structure", e.what());
            return;
        }catch(...){
            QMessageBox::critical(this, "Could not prepare structure", "Unrecognzied error when trying to prepare the structure.");
        }
    }else{
        vApp.newMol({vApp.curStep(), vApp.curMol->name + " (" + FFname + ')'});
    }
}

void LammpsWidget::mkGeom(const Step &curStep, const ForceField &FF, const fs::path &tempdir)
{
    // setup preset for FF demands
    auto preset = Plugins::LmpInput.makePreset();
    std::get<NamedEnum>(preset.at("style").first) = "full";
    preset.at("coeff").first = true;
    preset.at("bonds").first = FF.bond.has_value();
    preset.at("angles").first = FF.angle.has_value();
    preset.at("dihedrals").first = FF.dihedral.has_value();
    preset.at("impropers").first = FF.improper.has_value();
    // request parameter from FF
    auto param = FF.prepareParameters(curStep);
    // create input file
    writeFile((tempdir/"geom.lmp").string(), &Plugins::LmpInput, *vApp.curMol,
              master->curVP->moldata[vApp.curMol].curStep-1,
              param, preset);
}

void LammpsWidget::mkScript(const Step &curStep, const ForceField &FF, const fs::path &tempdir)
{
    std::ofstream script{tempdir/"input"};
    if(!script){
        throw IOError{"Could not write input script at "+tempdir.string()+"/input"};
    }
    // common block: setup FF and cell
    script << fmt::format(
        "units real\n"
        "atom_style full\n"
        "boundary {box} {box} {box}\n"
        "pair_style {pair}\n"
        "bond_style {bond}\n"
        "angle_style {angle}\n"
        "dihedral_style {dihedral}\n"
        "improper_style {improper}\n"
        "{extra}\n"
        "read_data {dir}/geom.lmp\n"
        "thermo_style custom step temp etotal ke pe ebond eangle edihed eimp evdwl ecoul\n",
    fmt::arg("box", curStep.hasCell() ? 'p' : 's'),
    fmt::arg("pair", FF.pair.has_value() ? FF.pair.value() : "none"),
    fmt::arg("bond", FF.bond.has_value() ? FF.bond.value() : "none"),
    fmt::arg("angle", FF.angle.has_value() ? FF.angle.value() : "none"),
    fmt::arg("dihedral", FF.dihedral.has_value() ? FF.dihedral.value() : "none"),
    fmt::arg("improper", FF.improper.has_value() ? FF.improper.value() : "none"),
    fmt::arg("extra", fmt::join(FF.extra_cmds, "\n")),
    fmt::arg("dir", tempdir.string()));
    // run
    if(ui->calcStack->currentIndex() == 0){
        // Minimization
        script << fmt::format(
            "thermo 1\n"
            "fix vipster all vipster 1\n"
            "min_style {}\n",
        ui->minSel->currentText().toStdString()
        );
    }else{
        // MD
        auto temp_target = ui->tempInput->text().toDouble();
        auto press_target = ui->pressInput->text().toDouble();
        auto dynstep = ui->stepInput->text().toUInt();
        auto equibstep = ui->equibInput->text().toUInt();
        auto reportstep = ui->reportInput->text().toUInt();
        auto timestep = ui->timeInput->text().toDouble();
        auto velocities = ui->velSel->currentIndex();
        auto rampstep = equibstep ? equibstep : dynstep/10;
        auto ensemble = ui->mdSel->currentText();
        std::map<QString, std::string> ensembleFmt{
            {"nve", "fix integrate all nve"},
            {"nvt", "fix integrate all nvt temp {0} {0} 100.0"},
            {"nph", "fix integrate all nph iso {1} {1} 100.0"},
            {"npt", "fix integrate all npt temp {0} {0} 100.0 iso {1} {1} 100.0"},
        };
        script << fmt::format(
            "timestep {timestep}\n"
            "#init velocities:\n{velocities}\n"
            "{fix_integrate}\n"
            "#equilibrate\n{equilibrate}\n"
            "thermo {freq}\n"
            "fix vipster all vipster {freq}\n",
        fmt::arg("timestep", timestep),
        // velocities: 0 -> create random, 1 -> heat ramp, else -> use existing
        fmt::arg("velocities", fmt::format(velocities == 0 ? "#random\nvelocity all create {0} 123465 mom yes rot yes" :
                                           velocities == 1 ? "#heat ramp\nfix equib all nvt temp 0.1 {0} 100.0\nrun {1}\nunfix equib" :
                                           "#use existing", temp_target, rampstep)),
        fmt::arg("fix_integrate", fmt::format(ensembleFmt.at(ensemble), temp_target, press_target)),
        fmt::arg("equilibrate", equibstep ? fmt::format("run {}", equibstep) : "#no"),
        fmt::arg("freq", reportstep)
        );
    }
}

void LammpsWidget::on_ffSel_currentIndexChanged(int index){
    if(index == ui->ffSel->count()-1){
        ui->customFrame->setVisible(true);
        ui->ffPrepare->setHidden(true);
    }else{
        ui->customFrame->setHidden(true);
        ui->ffPrepare->setVisible(true);
    }
}
