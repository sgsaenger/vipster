#include "scripthelp.h"
#include "ui_scripthelp.h"
#include "vipster/filter.h"

ScriptHelp::ScriptHelp(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ScriptHelp)
{
    ui->setupUi(this);
    ui->filterLabel->setText(Vipster::FilterAbout);
    setWindowTitle("Script grammar");
}

ScriptHelp::~ScriptHelp()
{
    delete ui;
}
