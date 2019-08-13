#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QColorDialog>
#include "periodictablewidget.h"
#include "ui_periodictablewidget.h"
#include "mainwindow.h"

using namespace Vipster;

template<typename T>
void PeriodicTableWidget::registerProperty(QWidget*, T Element::*)
{}

template<>
void PeriodicTableWidget::registerProperty(QWidget* w, float Element::* prop)
{
    connect(static_cast<QDoubleSpinBox*>(w),
            qOverload<double>(&QDoubleSpinBox::valueChanged), this,
            [prop, this](double newVal){
                currentEntry->*prop = static_cast<float>(newVal);
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
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
void PeriodicTableWidget::registerProperty(QWidget* w, unsigned int Element::* prop)
{
    connect(static_cast<QSpinBox*>(w),
            qOverload<int>(&QSpinBox::valueChanged), this,
            [prop, this](int newVal){
                currentEntry->*prop = static_cast<unsigned int>(newVal);
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
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
void PeriodicTableWidget::registerProperty(QWidget* w, std::string Element::* prop)
{
    connect(static_cast<QLineEdit*>(w),
            &QLineEdit::editingFinished, this,
            [prop, w, this](){
                currentEntry->*prop = static_cast<QLineEdit*>(w)->text().toStdString();
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
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
void PeriodicTableWidget::registerProperty(QWidget* w, ColVec Element::* prop)
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
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
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

PeriodicTableWidget::PeriodicTableWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PeriodicTableWidget)
{
    ui->setupUi(this);
    registerProperty(ui->mSel, &Element::m);
    registerProperty(ui->zSel, &Element::Z);
    registerProperty(ui->covSel, &Element::covr);
    registerProperty(ui->vdwSel, &Element::vdwr);
    registerProperty(ui->cpnlSel, &Element::CPNL);
    registerProperty(ui->cpppSel, &Element::CPPP);
    registerProperty(ui->pwppSel, &Element::PWPP);
    registerProperty(ui->colSel, &Element::col);
    registerProperty(ui->cutSel, &Element::bondcut);
    emit(currentEntryChanged());
}

PeriodicTableWidget::~PeriodicTableWidget()
{
    delete ui;
}

void PeriodicTableWidget::setEntry(QListWidgetItem *item)
{
    if(item && table){
        currentEntry = &table->at(item->text().toStdString());
    }else{
        currentEntry = nullptr;
    }
    emit(currentEntryChanged());
}

void PeriodicTableWidget::setTable(PeriodicTable* pse)
{
    this->table = pse;
    ui->pseList->clear();
    if(pse == &Vipster::pte){
        isGlobal = true;
    }
    if(pse){
        for(const auto& entry: *pse){
            ui->pseList->addItem(QString::fromStdString(entry.first));
        }
    }
}

void PeriodicTableWidget::updateWidget(guiChange_t change)
{
    if(!isGlobal){
        if((change & guiMolChanged) == guiMolChanged){
            setTable(master->curMol->pte.get());
        }else if(change & GuiChange::atoms){
            setTable(this->table);
        }
    }
}
