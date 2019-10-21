#include "definewidget.h"
#include "ui_definewidget.h"
#include "mainwindow.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QColorDialog>
#include <QInputDialog>

using namespace Vipster;

DefineWidget::DefineWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::DefineWidget)
{
    ui->setupUi(this);
    contextActions.push_back(new QAction{"Update group", ui->defTable});
    contextActions.back()->setDisabled(true);
    connect(contextActions.back(), &QAction::triggered, this, &DefineWidget::updateAction);
    contextActions.push_back(new QAction{"Delete group", ui->defTable});
    contextActions.back()->setDisabled(true);
    connect(contextActions.back(), &QAction::triggered, this, &DefineWidget::deleteAction);
    contextActions.push_back(new QAction{"Set as selection", ui->defTable});
    contextActions.back()->setDisabled(true);
    connect(contextActions.back(), &QAction::triggered, this, &DefineWidget::toSelAction);
    ui->defTable->addActions(contextActions);
}

DefineWidget::~DefineWidget()
{
    delete ui;
    for(auto* i: contextActions){
        delete i;
    }
}

void DefineWidget::updateWidget(Vipster::GUI::change_t change)
{
    if((change & GUI::stepChanged) == GUI::stepChanged){
        curStep = master->curStep;
        curState = &master->curVP->stepdata[curStep];
        defMap = &curState->def;
        fillTable();
    }else if(change & GUI::Change::definitions){
        fillTable();
    }else if(change & (GUI::Change::atoms|GUI::Change::settings)){
        for(auto& def: *defMap){
            def.second.second->update(&def.second.first,
                                      master->settings.atRadVdW.val,
                                      master->settings.atRadFac.val);
        }
    }
}

void DefineWidget::fillTable()
{
    auto& table = *ui->defTable;
    table.clearSelection();
    QSignalBlocker blockTable{ui->defTable};
    // setup table
    int i{0};
    table.clearContents();
    ui->defTable->setRowCount(defMap->size());
    for(auto& def: *defMap){
        // if contained in extras, this group is displayed
        table.setItem(i, 0, new QTableWidgetItem{});
        table.item(i, 0)->setFlags(Qt::ItemIsSelectable|
                                   Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        table.item(i, 0)->setCheckState(Qt::CheckState(
            (std::find(curState->extras.begin(), curState->extras.end(), def.second.second)
             != curState->extras.end())*2));
        // name
        table.setItem(i, 1, new QTableWidgetItem(QString::fromStdString(def.first)));
        // filter-str
        table.setItem(i, 2, new QTableWidgetItem(QString::fromStdString(
            def.second.first.getFilter().toStr())));
        // color button
        auto* but = new QPushButton("Select");
        const auto& color = def.second.second->color;
        but->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                           .arg(color[0]).arg(color[1]).arg(color[2]));
        connect(but, &QPushButton::clicked, this, &DefineWidget::colButton_clicked);
        table.setCellWidget(i, 3, but);
        i++;
    }
}

void DefineWidget::on_newButton_clicked()
{
    bool ok{false};
    auto filter = QInputDialog::getText(this, "Create new filtered group",
                                        "Enter filter for new group:",
                                        QLineEdit::Normal, QString(), &ok
                                       ).toStdString();
    if(!ok) return;
    try {
        auto sel = curStep->select(filter);
        auto name = QInputDialog::getText(this, "Create new filtered group",
                                          "Enter name for new group:",
                                          QLineEdit::Normal, QString(), &ok
                                         ).toStdString();
        if(!ok) return;
        auto [it, _] = defMap->insert_or_assign(name, std::pair{std::move(sel),
            std::make_shared<GUI::SelData>(master->globals)});
        it->second.second->update(&it->second.first,
            master->settings.atRadVdW.val, master->settings.atRadFac.val);
        it->second.second->color = defaultColors[defMap->size()%5];
        master->curVP->stepdata[curStep].extras.push_back(it->second.second);
        triggerUpdate(GUI::Change::definitions);
    }catch(const Error &e){
        QMessageBox msg{this};
        msg.setText(QString{e.what()});
        msg.exec();
    }catch(...){
        QMessageBox msg{this};
        msg.setText("Unknown error when parsing new filter");
        msg.exec();
    }
}

void DefineWidget::deleteAction()
{
    if(curSel == defMap->end()){
        throw Error{"DefineWidget: \"delete group\" triggered with invalid selection"};
    }
    auto posExtra = std::find(curState->extras.begin(), curState->extras.end(),
                              curSel->second.second);
    if(posExtra != curState->extras.end()) curState->extras.erase(posExtra);
    defMap->erase(curSel);
    triggerUpdate(GUI::Change::definitions);
}

void DefineWidget::on_fromSelButton_clicked()
{
    bool ok{false};
    auto tmp = QInputDialog::getText(this, "Copy selection to filtered group",
                                     "Enter name for new group:",
                                     QLineEdit::Normal, QString(), &ok
                                    ).toStdString();
    if(!ok) return;
    auto [it, _] = defMap->insert_or_assign(tmp, std::pair{*master->curSel,
        std::make_shared<GUI::SelData>(master->globals)});
    it->second.second->update(&it->second.first,
        master->settings.atRadVdW.val, master->settings.atRadFac.val);
    it->second.second->color = defaultColors[defMap->size()%5];
    master->curVP->stepdata[curStep].extras.push_back(it->second.second);
    triggerUpdate(GUI::Change::definitions);
}

void DefineWidget::toSelAction()
{
    if(curSel == defMap->end()){
        throw Error{"DefineWidget: \"to selection\" triggered with invalid selection"};
    }
    *master->curSel = curSel->second.first;
    triggerUpdate(GUI::Change::selection);
}

void DefineWidget::updateAction()
{
    if(curSel == defMap->end()){
        throw Error{"DefineWidget: \"update group\" triggered with invalid selection"};
    }
    auto& step = curSel->second.first;
    step.setFilter(step.getFilter());
    triggerUpdate(GUI::Change::definitions);
}

void DefineWidget::on_defTable_cellChanged(int row, int column)
{
    QTableWidgetItem *cell = ui->defTable->item(row, column);
    switch(column){
    case 0:
        // toggle visibility
        if(cell->checkState()){
            curState->extras.push_back(curSel->second.second);
        }else{
            auto pos = std::find(curState->extras.begin(),
                                 curState->extras.end(),
                                 curSel->second.second);
            curState->extras.erase(pos);
        }
        break;
    case 1:
        // change name
        if(defMap->find(cell->text().toStdString()) != defMap->end()){
            QMessageBox msg{this};
            msg.setText("Name \""+cell->text()+"\" is already in use.");
            msg.exec();
            QSignalBlocker block{ui->defTable};
            cell->setText(curSel->first.c_str());
            return;
        }else{
            auto node = defMap->extract(curSel);
            node.key() = cell->text().toStdString();
            defMap->insert(std::move(node));
        }
        break;
    case 2:
        // change filter
        try{
            auto filter = cell->text().toStdString();
            curSel->second.first.setFilter(filter);
            curSel->second.first.evaluateCache();
            curSel->second.second->update(&curSel->second.first,
                master->settings.atRadVdW.val, master->settings.atRadFac.val);
            triggerUpdate(GUI::Change::definitions);
        }catch(const Error &e){
            QMessageBox msg{this};
            msg.setText(QString{e.what()});
            msg.exec();
        }catch(...){
            QMessageBox msg{this};
            msg.setText("Unknown error when parsing new filter");
            msg.exec();
        }
        break;
    }
}

void DefineWidget::on_defTable_itemSelectionChanged()
{
    auto sel = ui->defTable->selectedItems();
    if(sel.empty()){
        curSel = defMap->end();
        for(auto* a: contextActions){
            a->setDisabled(true);
        }
    }else{
        const auto& curName = ui->defTable->item(sel[0]->row(), 1)->text();
        curSel = defMap->find(curName.toStdString());
        for(auto* a: contextActions){
            a->setEnabled(true);
        }
    }
}

void DefineWidget::colButton_clicked()
{
    auto& col = curSel->second.second->color;
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
    static_cast<QPushButton*>(sender())->setStyleSheet(
        QString("background-color: %1").arg(newCol.name()));
    triggerUpdate(GUI::Change::definitions);
}

void DefineWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString{"About filters"}, FilterAbout);
}
