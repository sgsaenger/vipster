#include "scripthelp.h"
#include "ui_scripthelp.h"

ScriptHelp::ScriptHelp(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ScriptHelp)
{
    ui->setupUi(this);
}

ScriptHelp::~ScriptHelp()
{
    delete ui;
}
