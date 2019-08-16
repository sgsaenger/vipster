#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QColorDialog>
#include <QMessageBox>
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
                currentElement->*prop = static_cast<float>(newVal);
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
            [prop, w, this](){
                QSignalBlocker block{w};
                if(currentElement){
                    w->setEnabled(true);
                    static_cast<QDoubleSpinBox*>(w)->setValue(
                                static_cast<double>(currentElement->*prop));
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
                currentElement->*prop = static_cast<unsigned int>(newVal);
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
            [prop, w, this](){
                QSignalBlocker block{w};
                if(currentElement){
                    w->setEnabled(true);
                    static_cast<QSpinBox*>(w)->setValue(
                                static_cast<int>(currentElement->*prop));
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
                currentElement->*prop = static_cast<QLineEdit*>(w)->text().toStdString();
                triggerUpdate(GuiChange::settings);
            }
    );
    connect(this, &PeriodicTableWidget::currentEntryChanged, this,
            [prop, w, this](){
                QSignalBlocker block{w};
                if(currentElement){
                    w->setEnabled(true);
                    static_cast<QLineEdit*>(w)->setText(
                                QString::fromStdString(currentElement->*prop));
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
                ColVec& col = currentElement->*prop;
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
                QSignalBlocker block{w};
                if(currentElement){
                    w->setEnabled(true);
                    const ColVec& col = currentElement->*prop;
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
        auto tmp = table->find(item->text().toStdString());
        if(tmp != table->end()){
            // if all objects are valid and we got a known type, set everything up
            currentName = &tmp->first;
            currentElement = &tmp->second;
            if(isGlobal){
                bool defElem = pte.find(*currentName) != pte.end();
                ui->defaultBut->setEnabled(defElem);
                ui->defaultBut->setVisible(defElem);
                ui->deleteBut->setEnabled(!defElem);
                ui->deleteBut->setVisible(!defElem);
            }else{
                ui->toGlobalBut->setEnabled(true);
                ui->fromGlobalBut->setEnabled(master->pte.find(*currentName)
                                              != master->pte.end());
            }
            emit(currentEntryChanged());
            return;
        }
    }
    // if something's not right, disable widget
    currentName = nullptr;
    currentElement = nullptr;
    if(isGlobal){
        ui->defaultBut->setDisabled(true);
        ui->deleteBut->setDisabled(true);
    }else{
        ui->fromGlobalBut->setDisabled(true);
        ui->toGlobalBut->setDisabled(true);
    }
    emit(currentEntryChanged());
}

void PeriodicTableWidget::setTable(PeriodicTable* pte)
{
    table = pte;
    ui->pteList->clear();
    if(pte == &master->pte){
        isGlobal = true;
    }else{
        isGlobal = false;
    }
    ui->toGlobalBut->setVisible(!isGlobal);
    ui->fromGlobalBut->setVisible(!isGlobal);
    ui->defaultBut->setVisible(isGlobal);
    ui->deleteBut->setVisible(isGlobal);
    if(pte){
        for(const auto& entry: *pte){
            ui->pteList->addItem(QString::fromStdString(entry.first));
        }
    }
}

void PeriodicTableWidget::changeEntry()
{
    if(sender() == ui->toGlobalBut){
        master->pte[*currentName] = *currentElement;
    }else if(sender() == ui->fromGlobalBut){
        *currentElement = master->pte.at(*currentName);
    }else if(sender() == ui->defaultBut){
        *currentElement = Vipster::pte.at(*currentName);
    }else if(sender() == ui->deleteBut){
        table->erase(table->find(*currentName));
    }else{
        throw Error{"PeriodicTable: changeEntry called with unknown sender"};
    }
    setEntry(ui->pteList->currentItem());
    triggerUpdate(GuiChange::settings);
}

void PeriodicTableWidget::updateWidget(guiChange_t change)
{
    if(isGlobal){
        // pull in changes to global table
        setTable(table);
    }else{
        if((change & guiMolChanged) == guiMolChanged){
            // load table of new molecule...
            setTable(master->curMol->pte.get());
        }else if(change & GuiChange::atoms){
            // ...or update current table
            setTable(this->table);
        }
    }
}

void PeriodicTableWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About periodic tables"), Vipster::PeriodicTableAbout);
}
