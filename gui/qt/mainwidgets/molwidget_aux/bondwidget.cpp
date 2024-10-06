#include <QMessageBox>
#include "bondwidget.h"
#include "ui_bondwidget.h"
#include "vipsterapplication.h"
#include "bonddelegate.h"

using namespace Vipster;

BondWidget::BondWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::BondWidget)
{
    ui->setupUi(this);

    // configure table
    ui->bondTable->setModel(&bondModel);
    ui->bondTable->setItemDelegateForColumn(3, new BondDelegate{});

    // Connect ui elements
    connect(ui->displayButton, &QPushButton::toggled, ui->frame, &QWidget::setVisible);
    connect(ui->bondHelpButton, &QPushButton::clicked, [&](){
        QMessageBox::information(this, QString("About bonds"), Vipster::BondsAbout);
    });
    connect(ui->bondModeBox, &QComboBox::currentIndexChanged, this, &BondWidget::bondModeHandler);
    connect(ui->bondSetButton, &QPushButton::pressed, this, &BondWidget::recalculateBonds);
    connect(ui->ovlpTable, &QTableWidget::itemSelectionChanged, this, &BondWidget::ovlpTableSelectionHandler);

    // Connect to app state changes
    connect(&vApp, &Application::activeStepChanged, this, &BondWidget::setActiveStep);
    connect(&vApp, &Application::stepChanged, this, &BondWidget::updateStep);

    // Hide UI elements until requested
    ui->frame->hide();
    ui->ovlpTable->hide();
    ui->ovlpLabel->hide();
}

BondWidget::~BondWidget()
{
    delete ui;
}

void BondWidget::setActiveStep(Step &step, Step::selection &)
{
    // expose BondMode
    QSignalBlocker blockBondMode(ui->bondModeBox);
    const bool autobonds = vApp.getState(step).automatic_bonds;
    ui->bondModeBox->setCurrentIndex(autobonds ? 1 : 0);
    ui->bondSetButton->setDisabled(autobonds);

    // enable bond widget when manual bonds are enabled
    if(!autobonds){
        ui->displayButton->setChecked(true);
    }

    // reset table model
    bondModel.reset();
}

void BondWidget::updateStep(Step &step)
{
    // only update active step
    if (&step != &vApp.curStep()) return;

    // reset table model
    bondModel.reset();

    updateOverlap();
}

void BondWidget::recalculateBonds()
{
    vApp.invokeOnStep(&Step::generateBonds, false);
}

void BondWidget::bondModeHandler(int index)
{
    auto is_automatic = static_cast<bool>(index);
    ui->bondSetButton->setDisabled(is_automatic);
    // TODO: this is a workaround to update app state, invoke on StepState instead of Step?
    vApp.invokeOnStep([](Step &s, bool is_automatic){
        vApp.getState(s).automatic_bonds = is_automatic;
    }, is_automatic);
}

void BondWidget::ovlpTableSelectionHandler()
{
    auto selection = ui->ovlpTable->selectedItems();
    if(!selection.empty()){
        const auto& sel = *selection[0];
        const auto& ovlp = vApp.curStep().getOverlaps()[sel.row()];
        const auto idx = sel.column() == 0 ? ovlp.at1 : ovlp.at2;
        vApp.invokeOnSel([](Step::selection &sel, size_t idx){
            SelectionFilter filter{};
            filter.mode = SelectionFilter::Mode::Index;
            filter.indices.emplace_back(idx, SizeVec{});
            // TODO: sort out const-correctness
            sel = const_cast<Step&>(vApp.curStep()).select(filter);
        }, idx);
    }
}

void BondWidget::updateOverlap()
{
    const auto& ovlp = vApp.curStep().getOverlaps();
    if(ovlp.empty()){
        ui->ovlpTable->hide();
        ui->ovlpLabel->hide();
        ui->displayButton->setText("Bonds");
    }else{
        auto& table = *ui->ovlpTable;
        QSignalBlocker tableblock{table};
        table.setCurrentCell(-1, -1);
        table.setRowCount(ovlp.size());
        for(size_t i=0; i<ovlp.size(); ++i){
            table.setItem(i, 0, new QTableWidgetItem{QString::number(ovlp[i].at1)});
            table.setItem(i, 1, new QTableWidgetItem{QString::number(ovlp[i].at2)});
        }
        table.show();
        ui->ovlpLabel->show();
        ui->displayButton->setText("Bonds (!)");
    }
}
