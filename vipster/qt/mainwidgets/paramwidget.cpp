#include "paramwidget.h"
#include "ui_paramwidget.h"

#include <QMessageBox>

using namespace Vipster;

ParamWidget::ParamWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ParamWidget)
{
    ui->setupUi(this);
    formats = makeParamWidgets();
    for(auto& p: formats){
        ui->paramStack->addWidget(p.second);
    }
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
                          "(" +  fmt->command +
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
    auto pos = formats.find(curFmt);
    if(pos == formats.end()){
        throw Error("Invalid parameter format");
    }
    curParam = pair.second.get();
    ui->paramStack->setCurrentWidget(pos->second);
    pos->second->setParam(curParam);
}

void ParamWidget::on_pushButton_clicked()
{
    QMessageBox::information(this, QString("About parameter presets"), Vipster::IO::ParametersAbout);
}
