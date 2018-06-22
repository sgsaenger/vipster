#include "paramwidget.h"
#include "ui_paramwidget.h"
#include "iowrapper.h"

using namespace Vipster;

ParamBase::ParamBase(QWidget *parent):
    QWidget(parent)
{}

ParamWidget::ParamWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ParamWidget)
{
    ui->setupUi(this);
}

ParamWidget::~ParamWidget()
{
    delete ui;
}

void ParamWidget::registerParam(Vipster::IOFmt fmt,
                                std::unique_ptr<Vipster::BaseParam>&& data)
{
    params.emplace_back(fmt, std::move(data));
    ui->paramSel->addItem(QString::fromStdString(
                          "(" +  IOPlugins.at(fmt)->command +
                           ") " + params.back().second->name
                         ));
}

void ParamWidget::on_paramSel_currentIndexChanged(int index)
{
    const auto& pair = params.at(static_cast<size_t>(index));
    curParam = pair.second.get();
    switch(pair.first){
    case IOFmt::PWI:
        ui->paramStack->setCurrentWidget(ui->PWWidget);
        ui->PWWidget->setParam(curParam);
        break;
    default:
        throw Error("Invalid parameter format");
    }
}
