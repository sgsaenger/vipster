#include "cpparam.h"
#include "ui_cpparam.h"

using namespace Vipster;
static const QString sections[]{
    "&INFO", "&CPMD", "&SYSTEM", "&PIMD",
    "&PATH", "&PTDDFT", "&ATOMS", "&DFT",
    "&PROP", "&RESP", "&LINRES", "&TDDFT",
    "&HARDNESS", "&CLASSIC", "&EXTE", "&VDW",
    "&QMMM",
};

CPParam::CPParam(QWidget *parent) :
    ParamBase(parent),
    ui(new Ui::CPParam)
{
    ui->setupUi(this);
    for(const auto& s: sections){
        ui->comboBox->addItem(s);
    }
}

CPParam::~CPParam()
{
    delete ui;
}

void CPParam::setParam(IO::Parameter *p)
{
    saveText();
    curParam = p;
    if(curParam->getFmt() != &IO::CPInput){
        throw Error("Invalid parameter set");
    }
    fillText();
    ui->prefixEdit->setText(std::get<std::string>(curParam->at("PPPrefix")).c_str());
    ui->suffixEdit->setText(std::get<std::string>(curParam->at("PPSuffix")).c_str());
    ui->nlEdit->setText(std::get<std::string>(curParam->at("PPNonlocality")).c_str());
}

void CPParam::on_comboBox_currentIndexChanged(const QString &arg)
{
    saveText();
    if(!curParam){
        return;
    }
    curSection = &std::get<std::vector<std::string>>(curParam->at(arg.toStdString()));
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
    for(const auto& line: *curSection){
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
        curSection->clear();
        for(const auto& line: tmp.split('\n')){
            curSection->push_back(line.toStdString());
        }
    }
}

void CPParam::on_prefixEdit_editingFinished()
{
    if(!curParam) return;
    std::get<std::string>(curParam->at("PPPrefix")) =
            ui->prefixEdit->text().toStdString();
}

void CPParam::on_suffixEdit_editingFinished()
{
    if(!curParam) return;
    std::get<std::string>(curParam->at("PPSuffix")) =
            ui->suffixEdit->text().toStdString();
}

void CPParam::on_nlEdit_editingFinished()
{
    if(!curParam) return;
    std::get<std::string>(curParam->at("PPNonlocality")) =
            ui->nlEdit->text().toStdString();
}
