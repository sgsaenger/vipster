#include "cpparam.h"
#include "ui_cpparam.h"

using namespace Vipster;

CPParam::CPParam(QWidget *parent) :
    ParamBase(parent),
    ui(new Ui::CPParam)
{
    ui->setupUi(this);
    for(const auto& i: IO::CPParam::str2section){
        ui->comboBox->addItem(QString::fromStdString(i.first));
    }
}

CPParam::~CPParam()
{
    delete ui;
}

void CPParam::setParam(IO::BaseParam *p)
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
    auto cmp = arg1.toStdString();
    curSection = std::find_if(IO::CPParam::str2section.begin(),
                              IO::CPParam::str2section.end(),
                              [&cmp](const decltype(IO::CPParam::str2section)::value_type& pair){
                                  return pair.first == cmp;
                              })->second;
    fillText();
}

void CPParam::on_plainTextEdit_textChanged()
{
    saveText();
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
    auto tmp = ui->plainTextEdit->toPlainText();
    if(tmp.size()){
        auto& section = curParam->*curSection;
        section.clear();
        for(const auto& line: tmp.split('\n')){
            section.push_back(line.toStdString());
        }
    }
}
