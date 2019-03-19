#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QColorDialog>
#include "psewidget.h"
#include "ui_psewidget.h"
#include "mainwindow.h"

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
    connect(this, &PSEWidget::currentEntryChanged, this,
            [prop, w, this](){
                if(currentEntry){
                    w->setEnabled(true);
                    QSignalBlocker block{w};
                    static_cast<QDoubleSpinBox*>(w)->setValue(
                                static_cast<double>(currentEntry->*prop));
                }else{
                    w->setDisabled(true);
                }
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
    connect(this, &PSEWidget::currentEntryChanged, this,
            [prop, w, this](){
                if(currentEntry){
                    w->setEnabled(true);
                    QSignalBlocker block{w};
                    static_cast<QSpinBox*>(w)->setValue(
                                static_cast<int>(currentEntry->*prop));
                }else{
                    w->setDisabled(true);
                }
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
    connect(this, &PSEWidget::currentEntryChanged, this,
            [prop, w, this](){
                if(currentEntry){
                    w->setEnabled(true);
                    QSignalBlocker block{w};
                    static_cast<QLineEdit*>(w)->setText(
                                QString::fromStdString(currentEntry->*prop));
                }else{
                    w->setDisabled(true);
                }
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
                if(!newCol.isValid()){
                    return;
                }
                col = {static_cast<uint8_t>(newCol.red()),
                       static_cast<uint8_t>(newCol.green()),
                       static_cast<uint8_t>(newCol.blue()),
                       static_cast<uint8_t>(newCol.alpha())};
                w->setStyleSheet(QString("background-color: %1").arg(newCol.name()));
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PSEWidget::currentEntryChanged, this,
            [prop, w, this](){
                if(currentEntry){
                    w->setEnabled(true);
                    const ColVec& col = currentEntry->*prop;
                    w->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                                            .arg(col[0]).arg(col[1]).arg(col[2]));
                }else{
                    w->setDisabled(true);
                }
            }
    );
}

PSEWidget::PSEWidget(QWidget *parent) :
    BaseWidget(parent),
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
    emit(currentEntryChanged());
}

PSEWidget::~PSEWidget()
{
    delete ui;
}

void PSEWidget::setEntry(QListWidgetItem *item)
{
    if(item && pse){
        currentEntry = &pse->at(item->text().toStdString());
    }else{
        currentEntry = nullptr;
    }
    emit(currentEntryChanged());
}

void PSEWidget::setPSE(PseMap* pse)
{
    this->pse = pse;
    ui->pseList->clear();
    if(pse == &Vipster::pse){
        isGlobal = true;
    }
    if(pse){
        for(const auto& entry: *pse){
            ui->pseList->addItem(QString::fromStdString(entry.first));
        }
    }
}

void PSEWidget::updateWidget(guiChange_t change)
{
    if(!isGlobal){
        if((change & guiMolChanged) == guiMolChanged){
            setPSE(master->curMol->pse.get());
        }else if(change & GuiChange::atoms){
            setPSE(this->pse);
        }
    }
}
