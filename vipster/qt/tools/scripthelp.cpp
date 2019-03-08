#include "scripthelp.h"
#include "ui_scripthelp.h"
#include "stepsel.h"

ScriptHelp::ScriptHelp(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ScriptHelp)
{
    ui->setupUi(this);
    ui->filterLabel->setText(Vipster::FilterAbout);
}

ScriptHelp::~ScriptHelp()
{
    delete ui;
}
