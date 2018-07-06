#include "cpparam.h"
#include "ui_cpparam.h"

using namespace Vipster;

CPParam::CPParam(QWidget *parent) :
    ParamBase(parent),
    ui(new Ui::CPParam)
{
    ui->setupUi(this);
    QSignalBlocker block{ui->comboBox};
    for(const auto& i: IO::CPParam::str2section){
        ui->comboBox->addItem(QString::fromStdString(i.first));
    }
}

CPParam::~CPParam()
{
    delete ui;
}

void CPParam::setParam(BaseParam *p)
{
    saveText();
    curParam = dynamic_cast<IO::CPParam*>(p);
    if(!curParam){
        throw Error("Invalid parameter set");
    }
    fillText();
}

void CPParam::on_comboBox_currentIndexChanged(const QString &arg1)
{
    saveText();
    curSection = IO::CPParam::str2section.at(arg1.toStdString());
    fillText();
}

void CPParam::fillText()
{
    if(!curParam || !curSection){
        return;
    }
    ui->plainTextEdit->clear();
    QStringList tmp{};
    for(const auto& line: curParam->*curSection){
        tmp.append(QString::fromStdString(line));
    }
    ui->plainTextEdit->setPlainText(tmp.join('\n'));
}

void CPParam::saveText()
{
    if(!curParam || !curSection){
        return;
    }
    auto tmp = ui->plainTextEdit->toPlainText().split('\n');
    auto& section = curParam->*curSection;
    section.clear();
    for(const auto& line: tmp){
        section.push_back(line.toStdString());
    }
}

void CPParam::focusOutEvent(QFocusEvent*)
{
    saveText();
}
