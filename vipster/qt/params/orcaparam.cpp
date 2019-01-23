#include "orcaparam.h"
#include "ui_orcaparam.h"

using namespace Vipster;

ORCAParam::ORCAParam(QWidget *parent) :
    ParamBase(parent),
    ui(new Ui::ORCAParam)
{
    ui->setupUi(this);
}

ORCAParam::~ORCAParam()
{
    delete ui;
}

void ORCAParam::setParam(IO::BaseParam *p)
{
    saveText();
    curParam = dynamic_cast<IO::OrcaParam*>(p);
    if(!curParam){
        throw Error("Invalid parameter set");
    }
    ui->plainTextEdit->clear();
    QStringList tmp{};
    for(const auto& line: curParam->header){
        tmp.append(QString::fromStdString(line));
    }
    ui->plainTextEdit->setPlainText(tmp.join('\n'));
}

void ORCAParam::on_plainTextEdit_textChanged()
{
    saveText();
}

void ORCAParam::saveText()
{
    if(!curParam){
        return;
    }
    auto tmp = ui->plainTextEdit->toPlainText();
    if(tmp.size()){
        auto& header = curParam->header;
        header.clear();
        for(const auto& line: tmp.split('\n')){
            header.push_back(line.toStdString());
        }
    }
}
