#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QColorDialog>
#include "psewidget.h"
#include "ui_psewidget.h"

using namespace Vipster;

template<typename T>
void PSEWidget::registerProperty(QWidget*, T PseEntry::*)
{}

template<>
void PSEWidget::registerProperty(QWidget* w, float PseEntry::* prop)
{
    connect(static_cast<QDoubleSpinBox*>(w),
            qOverload<double>(&QDoubleSpinBox::valueChanged), this,
            [prop, this](double newVal){
                currentEntry->*prop = static_cast<float>(newVal);
                triggerUpdate(GuiChange::settings);
            }
    );
}

template<>
void PSEWidget::registerProperty(QWidget* w, unsigned int PseEntry::* prop)
{
    connect(static_cast<QSpinBox*>(w),
            qOverload<int>(&QSpinBox::valueChanged), this,
            [prop, this](int newVal){
                currentEntry->*prop = static_cast<unsigned int>(newVal);
                triggerUpdate(GuiChange::settings);
            }
    );
}

template<>
void PSEWidget::registerProperty(QWidget* w, std::string PseEntry::* prop)
{
    connect(static_cast<QLineEdit*>(w),
            &QLineEdit::editingFinished, this,
            [prop, w, this](){
                currentEntry->*prop = static_cast<QLineEdit*>(w)->text().toStdString();
                triggerUpdate(GuiChange::settings);
            }
    );
}

template<>
void PSEWidget::registerProperty(QWidget* w, ColVec PseEntry::* prop)
{
    connect(static_cast<QPushButton*>(w),
            &QPushButton::clicked, this,
            [prop, w, this](){
                ColVec& col = currentEntry->*prop;
                auto oldCol = QColor::fromRgb(col[0], col[1], col[2], col[3]);
                auto newCol = QColorDialog::getColor(oldCol, this, QString{},
                                                     QColorDialog::ShowAlphaChannel);
                col = {static_cast<uint8_t>(newCol.red()),
                       static_cast<uint8_t>(newCol.green()),
                       static_cast<uint8_t>(newCol.blue()),
                       static_cast<uint8_t>(newCol.alpha())};
                w->setStyleSheet(QString("background-color: %1").arg(newCol.name()));
                triggerUpdate(GuiChange::settings);
            }
    );
}

PSEWidget::PSEWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PSEWidget)
{
    ui->setupUi(this);
    registerProperty(ui->mSel, &PseEntry::m);
    registerProperty(ui->zSel, &PseEntry::Z);
    registerProperty(ui->covSel, &PseEntry::covr);
    registerProperty(ui->vdwSel, &PseEntry::vdwr);
    registerProperty(ui->cpnlSel, &PseEntry::CPNL);
    registerProperty(ui->cpppSel, &PseEntry::CPPP);
    registerProperty(ui->pwppSel, &PseEntry::PWPP);
    registerProperty(ui->colSel, &PseEntry::col);
    registerProperty(ui->cutSel, &PseEntry::bondcut);
}

PSEWidget::~PSEWidget()
{
    delete ui;
}

void PSEWidget::setEntry(QListWidgetItem *item)
{
    QSignalBlocker blockAll{this};
    currentEntry = &pse.at(item->text().toStdString());
    const PseEntry& entry = *currentEntry;
    ui->zSel->setValue(static_cast<int>(entry.Z));
    ui->mSel->setValue(static_cast<double>(entry.m));
    ui->cutSel->setValue(static_cast<double>(entry.bondcut));
    ui->covSel->setValue(static_cast<double>(entry.covr));
    ui->vdwSel->setValue(static_cast<double>(entry.vdwr));
    ui->pwppSel->setText(QString::fromStdString(entry.PWPP));
    ui->cpppSel->setText(QString::fromStdString(entry.CPPP));
    ui->cpnlSel->setText(QString::fromStdString(entry.CPNL));
    const ColVec& col = entry.col;
    ui->colSel->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                              .arg(col[0], col[1], col[2]));
}
