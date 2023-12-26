#include "molwidget.h"
#include "ui_molwidget.h"
#include "molwidget_aux/newelement.h"
#include "vipsterapplication.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QMenu>

using namespace Vipster;


constexpr const char* inactiveKpoints[] = {"Gamma", "Monkhorst-Pack grid", "Discrete"};
constexpr const char* activeKpoints[] = {"Gamma (active)",
                                         "Monkhorst-Pack grid (active)",
                                         "Discrete (active)"};

MolWidget::MolWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);

    // Connect ui elements
    connect(ui->typeButton, &QPushButton::toggled, ui->typeContainer, &QWidget::setVisible);

    // TODO: hide all by default?
    ui->typeContainer->setVisible(false);
}

MolWidget::~MolWidget()
{
    delete ui;
}

void MolWidget::updateWidget(GUI::change_t change)
{
    if (change & GUI::Change::atoms) {
        ui->typeWidget->setTable(&vApp.curMol().getPTE());
    }
}

void MolWidget::on_clearTableButton_clicked()
{
    vApp.curMol().cleanPTE();
    ui->typeWidget->setTable(&vApp.curMol().getPTE());
}

void MolWidget::on_newElemButton_clicked()
{
    if(newelement(vApp.curMol().getPTE()).exec()){
        ui->typeWidget->setTable(&vApp.curMol().getPTE());
    }
}
