#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QColorDialog>
#include <QMessageBox>

#include "periodictablewidget.h"
#include "ui_periodictablewidget.h"

#include "periodictable_aux/newelement.h"

#include "vipsterapplication.h"

using namespace Vipster;

template<>
void PeriodicTableWidget::registerProperty(QWidget* w, double Element::* prop)
{
    connect(static_cast<QDoubleSpinBox*>(w),
            qOverload<double>(&QDoubleSpinBox::valueChanged), this,
            [prop, this](double newVal){
                if(isGlobal){
                    vApp.invokeOnConfig([prop, this](ConfigState &s, double newVal){
                        s.periodicTable.at(*currentName).*prop = newVal;
                    }, newVal);
                }else{
                    vApp.invokeOnStep([prop, this](Step &s, double newVal){
                        s.getPTE().at(*currentName).*prop = newVal;
                    }, newVal);
                }
            }
    );
    connect(this, &PeriodicTableWidget::currentElementChanged, w,
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
                if(isGlobal){
                    vApp.invokeOnConfig([prop, this](ConfigState &s, int newVal){
                        s.periodicTable.at(*currentName).*prop = newVal;
                    }, newVal);
                }else{
                    vApp.invokeOnStep([prop, this](Step &s, int newVal){
                        s.getPTE().at(*currentName).*prop = newVal;
                    }, newVal);
                }
            }
    );
    connect(this, &PeriodicTableWidget::currentElementChanged, this,
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
                auto newVal = static_cast<QLineEdit*>(w)->text().toStdString();
                if(isGlobal){
                    vApp.invokeOnConfig([prop, this](ConfigState &s, const std::string &newVal){
                        s.periodicTable.at(*currentName).*prop = newVal;
                    }, newVal);
                }else{
                    vApp.invokeOnStep([prop, this](Step &s, const std::string &newVal){
                        s.getPTE().at(*currentName).*prop = newVal;
                    }, newVal);
                }
            }
    );
    connect(this, &PeriodicTableWidget::currentElementChanged, this,
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
                const ColVec& col = currentElement->*prop;
                auto oldCol = QColor::fromRgb(col[0], col[1], col[2], col[3]);
                auto newCol = QColorDialog::getColor(oldCol, this, QString{},
                                                     QColorDialog::ShowAlphaChannel);
                if(!newCol.isValid()){
                    return;
                }
                ColVec newVal = {static_cast<uint8_t>(newCol.red()),
                                 static_cast<uint8_t>(newCol.green()),
                                 static_cast<uint8_t>(newCol.blue()),
                                 static_cast<uint8_t>(newCol.alpha())};
                if(isGlobal){
                    vApp.invokeOnConfig([prop, this](ConfigState &s, const ColVec &newVal){
                        s.periodicTable.at(*currentName).*prop = newVal;
                    }, newVal);
                }else{
                    vApp.invokeOnStep([prop, this](Step &s, const ColVec &newVal){
                        s.getPTE().at(*currentName).*prop = newVal;
                    }, newVal);
                }
                w->setStyleSheet(QString("background-color: %1").arg(newCol.name()));
            }
    );
    connect(this, &PeriodicTableWidget::currentElementChanged, this,
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

PeriodicTableWidget::PeriodicTableWidget(QWidget *parent, bool isGlob) :
    ui(new Ui::PeriodicTableWidget),
    isGlobal{isGlob}
{
    ui->setupUi(this);

    connect(ui->pteList, &QListWidget::currentItemChanged, this, &PeriodicTableWidget::setElement);

    // Connect ui elements for properties of PTE entries
    registerProperty(ui->mSel, &Element::m);
    registerProperty(ui->zSel, &Element::Z);
    registerProperty(ui->covSel, &Element::covr);
    registerProperty(ui->vdwSel, &Element::vdwr);
    registerProperty(ui->cpnlSel, &Element::CPNL);
    registerProperty(ui->cpppSel, &Element::CPPP);
    registerProperty(ui->pwppSel, &Element::PWPP);
    registerProperty(ui->colSel, &Element::col);
    registerProperty(ui->cutSel, &Element::bondcut);

    //initialize table
    if(isGlobal){
        connect(&vApp, &Application::configChanged, this, [this](const ConfigState &cfg){
            setTable(&cfg.periodicTable);
        });

        // ensure that all regular types are present, regardless of user settings
        vApp.invokeOnConfig([](ConfigState &c){
            for(const auto&[el, _]: Vipster::periodicTable){
                c.periodicTable.find_or_fallback(el);
            }
        });
    } else {
        connect(&vApp, &Application::activeMolChanged, this, [this](const Molecule &m){
            setTable(&m.getPTE());
        });
        connect(&vApp, &Application::activeStepChanged, this, [this](const Step &s){
            setTable(&s.getPTE());
        });
        connect(&vApp, &Application::molChanged, this, [this](const Molecule &m){
            if (&m == &vApp.curMol()) {
                setTable(&m.getPTE());
            }
        });
        connect(&vApp, &Application::stepChanged, this, [this](const Step &s){
            if (&s == &vApp.curStep()) {
                setTable(&s.getPTE());
            }
        });
    }

    // setup buttons

    // display help
    connect(ui->helpButton, &QPushButton::clicked, this, [this](){
        QMessageBox::information(this, QString("About periodic tables"), Vipster::PeriodicTableAbout);
    });

    // create a new atom type
    connect(ui->newElemButton, &QPushButton::clicked, this, [isGlob](){
        newelement(isGlob).exec();
    });

    // store custom element in global table
    ui->toGlobalButton->setVisible(!isGlobal);
    connect(ui->toGlobalButton, &QPushButton::clicked, this, [this](){
        vApp.invokeOnConfig([](ConfigState &c, const std::string &newName, const Element &newElement){
            auto &table = c.periodicTable;
            table[newName] = newElement;
        }, *currentName, *currentElement);
    });

    // restore custom element from global table
    ui->fromGlobalButton->setVisible(!isGlobal);
    connect(ui->fromGlobalButton, &QPushButton::clicked, this, [this](){
        vApp.invokeOnStep([](Step &s, const std::string &curName){
            auto &table = s.getPTE();
            table[curName] = vApp.config().periodicTable.at(curName);
        }, *currentName);
    });

    // remove unused elements from step's table
    ui->cleanButton->setVisible(!isGlobal);
    connect(ui->cleanButton, &QPushButton::clicked, this, [this](){
        vApp.invokeOnMol([](Molecule &m){
            m.cleanPTE();
        });
    });

    // restore global element to hard-coded default
    ui->defaultButton->setVisible(isGlobal);
    connect(ui->defaultButton, &QPushButton::clicked, this, [this](){
        vApp.invokeOnConfig([](ConfigState &c, const std::string &curName){
            c.periodicTable[curName] = Vipster::periodicTable.at(curName);
        }, *currentName);
    });

    // remove custom element
    ui->deleteButton->setVisible(isGlobal);
    connect(ui->deleteButton, &QPushButton::clicked, this, [this](){
        vApp.invokeOnConfig([](ConfigState &c, const std::string &curName){
            auto &table = c.periodicTable;
            table.erase(table.find(curName));
        }, *currentName);
    });
}

PeriodicTableWidget::~PeriodicTableWidget()
{
    delete ui;
}

void PeriodicTableWidget::setElement(QListWidgetItem *item)
{
    if(item && table){
        auto tmp = table->find(item->text().toStdString());
        if(tmp != table->end()){
            // if all objects are valid and we got a known type, set everything up
            currentName = &tmp->first;
            currentElement = &tmp->second;
            if(isGlobal){
                bool defElem = Vipster::periodicTable.find(*currentName) != Vipster::periodicTable.end();
                ui->defaultButton->setEnabled(defElem);
                ui->deleteButton->setEnabled(!defElem);
            }else{
                ui->toGlobalButton->setEnabled(true);
                const auto &globTable = vApp.config().periodicTable;
                ui->fromGlobalButton->setEnabled(globTable.find(*currentName)
                                              != globTable.end());
            }
            emit(currentElementChanged());
            return;
        }
    }
    // if no item is selected or something's not right, disable widget
    currentName = nullptr;
    currentElement = nullptr;
    if(isGlobal){
        ui->defaultButton->setDisabled(true);
        ui->deleteButton->setDisabled(true);
    }else{
        ui->fromGlobalButton->setDisabled(true);
        ui->toGlobalButton->setDisabled(true);
    }
    emit(currentElementChanged());
}

void PeriodicTableWidget::setTable(const PeriodicTable* pte)
{
    table = pte;
    auto row = ui->pteList->currentRow();
    ui->pteList->clear();
    if(pte){
        for(const auto& entry: *pte){
            ui->pteList->addItem(QString::fromStdString(entry.first));
        }
    }
    if (row > ui->pteList->count()) row = ui->pteList->count()-1;
    ui->pteList->setCurrentRow(row);
}
