#include "lammpswidget.lmp.h"
#include "ui_lammpswidget.lmp.h"
#include "../mainwindow.h"
#include "io/plugins/lmpinput.h"

#include <fmt/format.h>

#include <QMessageBox>

#include "lammps/lammps.h"
#include "lammps/atom.h"
#include "lammps/comm.h"
#include "lammps/domain.h"
#include "lammps/error.h"
#include "lammps/input.h"
#include "lammps/info.h"
#include "lammps/update.h"
#include "lammps/exceptions.h"

using namespace Vipster;
using namespace LAMMPS_NS;
namespace fs = std::filesystem;

LammpsWidget::LammpsWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::LammpsWidget),
    forcefields{defaultForcefields()}
{
    ui->setupUi(this);
    ui->customFrame->setHidden(true);
    // set validators for minimize settings
    ui->etolInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->ftolInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->maxItInput->setValidator(new QIntValidator{1, 10000000});
    ui->maxevalInput->setValidator(new QIntValidator{1, 10000000});
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
    LAMMPS lmp{0, nullptr, MPI_COMM_WORLD};
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

    // TODO: check if we can register callback?
}

LammpsWidget::~LammpsWidget()
{
    delete ui;
}

fs::path getLmpTmpDir()
{
    return getTempPath()/"lammps";
}

void LammpsWidget::on_runButton_clicked()
{
    if(!master->curStep->getNat()){
        return;
    }
    // get forcefield
    const auto FFname = ui->ffSel->currentText().toStdString();
    if(FFname == "Custom"){
        return;
    }
    const auto& FF = defaultForcefields().at(FFname);
    // launch lammps in tempdir
    auto tempdir = getLmpTmpDir();
    auto log = tempdir/"log.lammps";
    fs::create_directory(tempdir);
    std::array<char*, 5> lmparg{
        nullptr,
        "-screen", "none",
        "-log", &log.string()[0]
    };
    LAMMPS lmp{5, lmparg.data(), MPI_COMM_WORLD};
    try{
        auto preset = IO::LmpInput.makePreset();
        // TODO: get this from forcefield?
//        preset.at("style").first = "full";
        std::get<bool>(preset.at("bonds").first) = FF->required_bond.has_value();
//        preset.at("angles").first = FF->required_angle.has_value();
//        preset.at("dihedrals").first = FF->required_dihedral.has_value();
//        preset.at("impropers").first = FF->required_improper.has_value();
        // FIXME: LmpInput is BROKEN. maybe because selection shares bonds with step which it SHOULD NOT FIX IT FIX IT
        writeFile(tempdir/"inpgeom", &IO::LmpInput, *master->curMol,
                  master->curVP->moldata[master->curMol].curStep-1,
                  std::nullopt, preset);
        auto curStep = master->curStep->asFmt(AtomFmt::Angstrom);
        // initialize lammps
        lmp.init();
        // always use real units -> angstrom
        lmp.update->set_units("real");
        // always use full atoms because why not
        lmp.input->one("atom_style full");
        // create box
        Mat vec;
        if(curStep.hasCell()){
            // use existing, periodic cell
            lmp.input->one("boundary p p p");
            vec = curStep.getCellVec();
            if(!float_comp(vec[0][1], 0.) || !float_comp(vec[0][2], 0.) || !float_comp(vec[1][2], 0.)){
                QMessageBox::critical(this, "Error in LAMMPS setup",
                                      "LAMMPS: Cell vectors must form diagonal or lower triangular matrix");
                return;
            }
            if(!float_comp(vec[1][0], 0.) || !float_comp(vec[2][0], 0.) || !float_comp(vec[2][1], 0.)){
                lmp.input->one(fmt::format("region box prism {} {} {} {} {} {} {} {} {}",
                                           0, vec[0][0],
                                           0, vec[1][1],
                                           0, vec[2][2],
                                           vec[1][0], vec[2][0], vec[2][1]).c_str());
            }else{
                lmp.input->one(fmt::format("region box block {} {} {} {} {} {}",
                                           0, vec[0][0],
                                           0, vec[1][1],
                                           0, vec[2][2]).c_str());
            }
        }else{
            // create a shrink-wrapped box for lammps
            lmp.input->one("boundary s s s");
            Vec pos_min{{std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()}};
            Vec pos_max{{std::numeric_limits<double>::lowest(),
                     std::numeric_limits<double>::lowest(),
                     std::numeric_limits<double>::lowest()}};
            for(const auto& at: curStep){
                pos_min[0] = std::min(pos_min[0], at.coord[0]);
                pos_min[1] = std::min(pos_min[1], at.coord[1]);
                pos_min[2] = std::min(pos_min[2], at.coord[2]);
                pos_max[0] = std::max(pos_max[0], at.coord[0]);
                pos_max[1] = std::max(pos_max[1], at.coord[1]);
                pos_max[2] = std::max(pos_max[2], at.coord[2]);
            }
            lmp.input->one(fmt::format("region box block {} {} {} {} {} {}",
                                       pos_min[0], pos_max[0],
                                       pos_min[1], pos_max[1],
                                       pos_min[2], pos_max[2]).c_str());
        }
        lmp.input->one(fmt::format("create_box {} box", curStep.getNtyp()).c_str());
    }catch(LAMMPSAbortException &e){
        QMessageBox::critical(this, "Error in LAMMPS run", QString::fromStdString(e.message));
    }catch(LAMMPSException &e){
        QMessageBox::warning(this, "Error in LAMMPS run", QString::fromStdString(e.message));
    }
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
    const auto& FF = defaultForcefields().at(FFname);
    if(FF->prepareStep){
        try{
            auto mol = FF->prepareStep(*master->curStep, master->curMol->name);
            master->newMol(std::move(mol));
        }catch(const Vipster::Error &e){
            QMessageBox::critical(this, "Could not deduce atom type", e.what());
            return;
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
