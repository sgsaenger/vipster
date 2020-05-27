#include "newelement.h"
#include "ui_newelement.h"
#include <QMessageBox>

using namespace Vipster;

newelement::newelement(Vipster::PeriodicTable &table, QWidget *parent) :
    QDialog{parent},
    ui{new Ui::newelement},
    table{table}
{
    ui->setupUi(this);
    setWindowTitle("Create or replace Element");
}

newelement::~newelement()
{
    delete ui;
}

void newelement::accept()
{
    auto newName = ui->newName->text().toStdString();
    if(ui->baseCheck->isChecked()){
        auto baseName = ui->baseName->text().toStdString();
        auto baseElem = table.root->find(baseName);
        if(baseElem == table.root->end()){
            QMessageBox::critical(this, "Invalid base element", QString{"The element \""}+baseName.c_str()+
                                  "\" is not known. Please choose a valid base element.");
            return;
        }
        table[newName] = baseElem->second;
    }else{
        table.find_or_fallback(newName);
    }
    QDialog::accept();
}
