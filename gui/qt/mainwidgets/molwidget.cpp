#include "molwidget.h"
#include "ui_molwidget.h"
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

    // Periodictable Widget is used also for global table
    // -> containerization handled here
    connect(ui->typeButton, &QPushButton::toggled, ui->typeContainer, &QWidget::setVisible);
    ui->typeContainer->setVisible(ui->typeButton->isChecked());
}

MolWidget::~MolWidget()
{
    delete ui;
}
