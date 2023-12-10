#include <QMessageBox>
#include "cellwidget.h"
#include "vipsterapplication.h"
#include "doubledelegate.h"
#include "ui_cellwidget.h"

using namespace Vipster;

CellWidget::CellWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::CellWidget)
{
    ui->setupUi(this);

    // setup cell-table
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
             ui->cellVecTable->setItem(j,k,new QTableWidgetItem());
        }
    }
    ui->cellVecTable->setItemDelegate(new DoubleDelegate{});

    // Connect ui elements
    connect(ui->displayButton, &QPushButton::toggled, ui->frame, &QWidget::setVisible);
    connect(ui->cellEnabledBox, &QCheckBox::toggled, this, &CellWidget::enableCell);
    connect(ui->cellFmt, &QComboBox::currentIndexChanged, this, &CellWidget::cellFmtChanged);
    connect(ui->cellDimBox, &QDoubleSpinBox::valueChanged, this, &CellWidget::cellDimChanged);
    connect(ui->cellVecTable, &QTableWidget::cellChanged, this, &CellWidget::cellVecChanged);
    connect(ui->cellTrajecButton, &QPushButton::clicked, this, &CellWidget::applyToTrajectory);

    // Connect to app state changes
    connect(&vApp, &Application::activeStepChanged, this, &CellWidget::updateStep);
    connect(&vApp, &Application::stepChanged, this, &CellWidget::updateStep);

    // hide ui until requested
    ui->frame->hide();
}

CellWidget::~CellWidget()
{
    delete ui;
}

void CellWidget::updateStep(Step &step)
{
    // only update active step
    if (&step != &vApp.curStep()) return;

    QSignalBlocker blockEnabled(ui->cellEnabledBox);
    ui->cellEnabledBox->setChecked(step.hasCell());
    if(step.hasCell()){
        ui->displayButton->setChecked(true);
    }

    QSignalBlocker blockDim(ui->cellDimBox);
    ui->cellDimBox->setValue(static_cast<double>(
        step.getCellDim(static_cast<AtomFmt>(ui->cellFmt->currentIndex()))
    ));

    QSignalBlocker blockTable(ui->cellVecTable);
    const Mat& vec = step.getCellVec();
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
            ui->cellVecTable->item(j,k)->setText(QString::number(vec[j][k]));
        }
    }
}

void CellWidget::enableCell(bool enabled)
{
    if(enabled){
        const auto scale = ui->cellScaleBox->isChecked();
        const auto dim = ui->cellDimBox->value();
        const auto fmt = static_cast<AtomFmt>(ui->cellFmt->currentIndex());
        const auto& table = *ui->cellVecTable;
        Mat cell{
            Vec{
                table.item(0,0)->text().toDouble(),
                table.item(0,1)->text().toDouble(),
                table.item(0,2)->text().toDouble(),
            }, Vec{
                table.item(1,0)->text().toDouble(),
                table.item(1,1)->text().toDouble(),
                table.item(1,2)->text().toDouble(),
            }, Vec{
                table.item(2,0)->text().toDouble(),
                table.item(2,1)->text().toDouble(),
                table.item(2,2)->text().toDouble(),
            }
        };
        try{
            if(cell == Mat{}){
                vApp.invokeOnStep(
                    [](Step &step, const double &dim, AtomFmt fmt, bool scale){
                        step.enableCell(true);
                        step.setCellDim(dim, fmt, scale);
                    }, dim, fmt, scale);
            }else{
                vApp.invokeOnStep(
                    [](Step &step, const Mat &cell, const double &dim, AtomFmt fmt, bool scale){
                        step.setCellVec(cell, scale);
                        step.setCellDim(dim, fmt, scale);
                    }, cell, dim, fmt, scale);
            }
        } catch(const Error& e){
            QMessageBox::critical(this, "Error setting cell vectors", e.what());
            QSignalBlocker block{ui->cellEnabledBox};
            ui->cellEnabledBox->setCheckState(Qt::CheckState::Unchecked);
        }
    }else{
        vApp.invokeOnStep(&Step::enableCell, false);
    }
}

void CellWidget::cellFmtChanged(int idx)
{
    QSignalBlocker blockCDB(ui->cellDimBox);
    ui->cellDimBox->setValue(static_cast<double>(vApp.curStep().getCellDim(static_cast<AtomFmt>(idx))));
}

void CellWidget::cellDimChanged(double dim)
{
    // if cell is disabled, exit early
    if(!ui->cellEnabledBox->isChecked()){
        return;
    }

    const auto fmt = static_cast<AtomFmt>(ui->cellFmt->currentIndex());
    const auto scale = ui->cellScaleBox->isChecked();
    vApp.invokeOnStep(&Step::setCellDim, dim, fmt, scale);
}

void CellWidget::cellVecChanged(int row, int col)
{
    // if cell is disabled, exit early
    if(!ui->cellEnabledBox->isChecked()){
        return;
    }

    const auto scale = ui->cellScaleBox->isChecked();
    const auto val = ui->cellVecTable->item(row, col)->text().toDouble();
    Mat cell = vApp.curStep().getCellVec();
    const auto origVal = cell[row][col];
    cell[row][col] = val;
    try {
        vApp.invokeOnStep(&Step::setCellVec, cell, scale);
    } catch (const Error &e) {
        QMessageBox::critical(this, "Error setting cell vectors", e.what());
        QSignalBlocker block{ui->cellVecTable};
        ui->cellVecTable->item(row, col)->setText(QString::number(origVal));
    }
}

void CellWidget::applyToTrajectory()
{
    if(ui->cellEnabledBox->isChecked()){
        const auto scale = ui->cellScaleBox->isChecked();
        const auto dim = ui->cellDimBox->value();
        const auto fmt = static_cast<AtomFmt>(ui->cellFmt->currentIndex());
        const auto cell = vApp.curStep().getCellVec();
        vApp.invokeOnTrajec(
            [](Step &step, const Mat &cell, const double &dim, AtomFmt fmt, bool scale){
                step.setCellVec(cell, scale);
                step.setCellDim(dim, fmt, scale);
            }, cell, dim, fmt, scale);
    }else{
        vApp.invokeOnTrajec(&Step::enableCell, false);
    }
}
