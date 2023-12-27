#include "newelement.h"
#include "ui_newelement.h"
#include <QMessageBox>

#include "vipsterapplication.h"

#include "vipster/molecule.h"

using namespace Vipster;

newelement::newelement(bool isGlobal, QWidget *parent) :
    QDialog{parent},
    ui{new Ui::newelement},
    isGlobal{isGlobal}
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
    auto baseName = ui->baseCheck->isChecked()
                    ? std::make_optional<std::string>(ui->baseName->text().toStdString())
                    : std::nullopt;
    if (baseName) {
        try {
            if (isGlobal) {
                vApp.invokeOnConfig([&](ConfigState &c){
                    auto &table = c.periodicTable;
                    table[newName] = table.at(*baseName);
                });
            } else {
                vApp.invokeOnStep([&](Step &s){
                    auto &table = s.getPTE();
                    table[newName] = table.root->at(*baseName);
                });
            }
        } catch (const std::out_of_range &) {
            QMessageBox::critical(this, "Invalid base element",
                                  QString{"The element \""} + baseName->c_str() +
                                  "\" is not known. Please choose a valid base element.");
            return;
        }
    } else {
        if (isGlobal) {
            vApp.invokeOnConfig([&](ConfigState &c){
                c.periodicTable.find_or_fallback(newName);
            });
        } else {
            vApp.invokeOnStep([&](Step &s){
                s.getPTE().find_or_fallback(newName);
            });
        }
    }
    QDialog::accept();
}
