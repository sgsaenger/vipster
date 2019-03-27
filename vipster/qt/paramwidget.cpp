#include "paramwidget.h"
#include "ui_paramwidget.h"

#include <QMessageBox>

using namespace Vipster;

ParamBase::ParamBase(QWidget *parent):
    QWidget(parent)
{}

ParamWidget::ParamWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ParamWidget)
{
    ui->setupUi(this);
}

ParamWidget::~ParamWidget()
{
    delete ui;
}

void ParamWidget::clearParams()
{
    params.clear();
    ui->paramSel->clear();
}

void ParamWidget::registerParam(std::unique_ptr<Vipster::IO::BaseParam>&& data)
{
    auto fmt = data->getFmt();
    params.emplace_back(fmt, std::move(data));
    ui->paramSel->addItem(QString::fromStdString(
                          "(" +  IOPlugins.at(fmt)->command +
                           ") " + params.back().second->name
                         ));
    ui->paramSel->setCurrentIndex(ui->paramSel->count()-1);
}

void ParamWidget::on_paramSel_currentIndexChanged(int index)
{
    if(index<0){
        ui->paramStack->setCurrentWidget(ui->NoPWidget);
        curParam = nullptr;
        return;
    }
    if(static_cast<size_t>(index) >= params.size()){
        throw Error("Invalid parameter set selected");
    }
    const auto& pair = params.at(static_cast<size_t>(index));
    curFmt = pair.first;
    curParam = pair.second.get();
    switch(curFmt){
    case IOFmt::PWI:
        ui->paramStack->setCurrentWidget(ui->PWWidget);
        ui->PWWidget->setParam(curParam);
        break;
    case IOFmt::CPI:
        ui->paramStack->setCurrentWidget(ui->CPWidget);
        ui->CPWidget->setParam(curParam);
        break;
    case IOFmt::ORCA:
        ui->paramStack->setCurrentWidget(ui->ORCAWidget);
        ui->ORCAWidget->setParam(curParam);
        break;
    default:
        throw Error("Invalid parameter format");
    }
}

void ParamWidget::on_pushButton_clicked()
{
    QMessageBox::information(this, QString("About parameter presets"), Vipster::IO::ParametersAbout);
}
