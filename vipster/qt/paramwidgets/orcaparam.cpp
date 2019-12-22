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

void ORCAParam::setParam(IO::Parameter *p)
{
    saveText();
    curParam = p;
    if(curParam->getFmt() != &IO::OrcaInput){
        throw Error("Invalid parameter set");
    }
    ui->plainTextEdit->clear();
    QStringList tmp{};
    for(const auto& line: std::get<std::vector<std::string>>(curParam->at("header").first)){
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
        auto& header = std::get<std::vector<std::string>>(curParam->at("header").first);
        header.clear();
        for(const auto& line: tmp.split('\n')){
            header.push_back(line.toStdString());
        }
    }
}
