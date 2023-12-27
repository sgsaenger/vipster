#include "molwidget.h"
#include "ui_molwidget.h"
#include "vipsterapplication.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QMenu>

using namespace Vipster;

MolWidget::MolWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);

    // Connect ui elements
    connect(ui->typeButton, &QPushButton::toggled, ui->typeContainer, &QWidget::setVisible);
}

MolWidget::~MolWidget()
{
    delete ui;
}
