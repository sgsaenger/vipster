#include "molwidget.h"
#include "ui_molwidget.h"

using namespace Vipster;

MolWidget::MolWidget(QWidget *parent) :
    QScrollArea(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);

    // Periodictable Widget is used also for global table
    // -> containerization handled here
    connect(ui->typeButton, &QPushButton::toggled, ui->typeContainer, &QWidget::setVisible);
    ui->typeContainer->setVisible(ui->typeButton->isChecked());
}

MolWidget::~MolWidget()
{
    delete ui;
}
