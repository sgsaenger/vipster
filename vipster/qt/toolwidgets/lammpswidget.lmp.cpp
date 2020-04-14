#include "lammpswidget.lmp.h"
#include "ui_lammpswidget.lmp.h"
#include "../mainwindow.h"
#include "io/plugins/lmpinput.h"
#include "lammpswidget_aux/fix_vipster.lmp.h"

#include <fmt/format.h>

#include <QMessageBox>

#include "lammps/lammps.h"
#include "lammps/atom.h"
#include "lammps/comm.h"
#include "lammps/domain.h"
#include "lammps/error.h"
#include "lammps/input.h"
#include "lammps/info.h"
#include "lammps/modify.h"
#include "lammps/update.h"
#include "lammps/exceptions.h"

using namespace Vipster;
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
#ifdef USE_MPI
    ui->MPISpin->setDisabled(true);
#else
    ui->MPISpin->setEnabled(true);
#endif
    ui->OMPSpin->setDisabled(true);
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
        MPI_Init(&argc, &argv);
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
    delete ui;
}

fs::path getLmpTmpDir()
{
    // TODO: add random token
    return getTempPath()/"lammps";
}

void LammpsWidget::on_runButton_clicked()
{
    if(!master->curStep->getNat()){
        return;
    }
    auto curStep = master->curStep->asFmt(AtomFmt::Angstrom);
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
        // launch lammps in tempdir
        auto tempdir = getLmpTmpDir();
        auto log = (tempdir/"log.lammps").string();
        fs::create_directory(tempdir);
        std::array<char*, 5> lmparg{
            nullptr,
            "-screen", "none",
            "-log", &log[0]
        };
        LAMMPS lmp{5, lmparg.data(), MPI_COMM_WORLD};
        // setup preset for FF demands
        auto preset = IO::LmpInput.makePreset();
        std::get<NamedEnum>(preset.at("style").first) = "full";
        preset.at("coeff").first = true;
        preset.at("bonds").first = FF->bond.has_value();
        preset.at("angles").first = FF->angle.has_value();
        preset.at("dihedrals").first = FF->dihedral.has_value();
        preset.at("impropers").first = FF->improper.has_value();
        // request parameter from FF
        auto param = FF->prepareParameters(curStep);
        // create input file
        writeFile(tempdir/"inpgeom", &IO::LmpInput, *master->curMol,
                  master->curVP->moldata[master->curMol].curStep-1,
                  param, preset);
        // setup system for forcefield calculation
        lmp.input->one("units real"); // always use real units -> angstrom
        lmp.input->one("atom_style full"); // always use full atoms
        if(curStep.hasCell()){
            lmp.input->one("boundary p p p"); // periodic when cell is defined
        }else{
            lmp.input->one("boundary s s s"); // shrink-wrapped should create less surprises than fixed
        }
        if(FF->pair.has_value()){
            lmp.input->one(fmt::format("pair_style {}", FF->pair.value()).c_str());
        }
        if(FF->bond.has_value()){
            lmp.input->one(fmt::format("bond_style {}", FF->bond.value()).c_str());
        }
        if(FF->angle.has_value()){
            lmp.input->one(fmt::format("angle_style {}", FF->angle.value()).c_str());
        }
        if(FF->dihedral.has_value()){
            lmp.input->one(fmt::format("dihedral_style {}", FF->dihedral.value()).c_str());
        }
        if(FF->improper.has_value()){
            lmp.input->one(fmt::format("improper_style {}", FF->improper.value()).c_str());
        }
        for(const auto &cmd: FF->extra_cmds){
            lmp.input->one(cmd.c_str());
        }
        // read input file
        lmp.input->one(fmt::format("read_data {}/inpgeom", tempdir.string()).c_str());
        // get more details about energies
        lmp.input->one("thermo_style custom step temp etotal ke pe ebond eangle edihed eimp evdwl ecoul");
        // register custom fix
        (*lmp.modify->fix_map)["vipster"] = &mkFixVipster;
        // run
        if(ui->calcStack->currentIndex() == 0){
            // Minimization
            auto min_style = ui->minSel->currentText().toStdString();
            // report back each step
            lmp.input->one("thermo 1");
            lmp.input->one("fix vipster all vipster 1");
            auto fix_vipster = dynamic_cast<FixVipster*>(lmp.modify->fix[lmp.modify->nfix-1]);
            if(!fix_vipster){
                throw Vipster::Error{"Error on registering callback fix."};
            }else{
                fix_vipster->init_vipster(master, fmt::format("(Min: {})", min_style));
            }
            // set minimizer
            lmp.input->one(fmt::format("min_style {}", min_style).c_str());
            // trigger calculation
            lmp.input->one(fmt::format("minimize {} {} {} {}",
                                       ui->etolInput->text().toInt(),
                                       ui->ftolInput->text().toInt(),
                                       ui->maxItInput->text().toInt(),
                                       ui->maxevalInput->text().toInt()
                                       ).c_str());
        }else{
            // MD
            auto temp_target = ui->tempInput->text().toDouble();
            auto press_target = ui->pressInput->text().toDouble();
            auto dynstep = ui->stepInput->text().toUInt();
            auto equibstep = ui->equibInput->text().toUInt();
            auto reportstep = ui->reportInput->text().toUInt();
            // set timestep
            auto timestep = ui->timeInput->text().toDouble();
            lmp.input->one(fmt::format("timestep {}", timestep).c_str());
            // initialize velocities
            auto velocities = ui->velSel->currentIndex();
            if(velocities == 0){ // init random
                lmp.input->one(fmt::format("velocity all create {} 123465 mom yes rot yes", temp_target).c_str());
            }else if(velocities == 1){ // heat ramp
                auto rampstep = equibstep ? equibstep : dynstep/10;
                lmp.input->one(fmt::format("fix equib all nvt temp 0.1 {} 100.0", temp_target).c_str());
                lmp.input->one(fmt::format("run {}", rampstep).c_str());
                lmp.input->one("unfix equib");
            }// else use existing -> no preparation
            // integrating fix for MD
            auto ensemble = ui->mdSel->currentText();
            std::map<QString, std::string> ensembleFmt{
                {"nve", "fix integrate all nve"},
                {"nvt", "fix integrate all nvt temp {0} {0} 100.0"},
                {"nph", "fix integrate all nph iso {1} {1} 100.0"},
                {"npt", "fix integrate all npt temp {0} {0} 100.0 iso {1} {1} 100.0"},
            };
            lmp.input->one(fmt::format(ensembleFmt.at(ensemble), temp_target, press_target).c_str());
            // equilibration
            if(equibstep){
                lmp.input->one(fmt::format("run {}", equibstep).c_str());
            }
            // actual MD
            lmp.input->one(fmt::format("thermo {}", reportstep).c_str());
            lmp.input->one(fmt::format("fix vipster all vipster {}", reportstep).c_str());
            auto fix_vipster = dynamic_cast<FixVipster*>(lmp.modify->fix[lmp.modify->nfix-1]);
            if(!fix_vipster){
                throw Vipster::Error{"Error on registering callback fix."};
            }else{
                fix_vipster->init_vipster(master, fmt::format("(MD: {})", ensemble.toUpper().toStdString()));
            }
            lmp.input->one(fmt::format("run {}", dynstep).c_str());
        }
    }catch(LAMMPSAbortException &e){
        QMessageBox::critical(this, "Error in LAMMPS run", QString::fromStdString(e.what()));
    }catch(LAMMPSException &e){
        QMessageBox::warning(this, "Error in LAMMPS run", QString::fromStdString(e.what()));
    }catch(std::exception &e){
        QMessageBox::warning(this, "Error in LAMMPS run", QString::fromStdString(e.what()));
    }catch(...){
        QMessageBox::critical(this, "Error in LAMMPS run", "Unrecognized error when trying to run LAMMPS");
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
            auto mol = FF->prepareStep(*master->curStep, master->curMol->name);
            master->newMol(std::move(mol));
        }catch(const Vipster::Error &e){
            QMessageBox::warning(this, "Could not prepare structure", e.what());
            return;
        }catch(...){
            QMessageBox::critical(this, "Could not prepare structure", "Unrecognzied error when trying to prepare the structure.");
        }
    }else{
        master->newMol({*master->curStep, master->curMol->name + " (" + FFname + ')'});
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
