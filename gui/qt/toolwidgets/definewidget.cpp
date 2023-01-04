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
        curStep = vApp.curStep;
        defMap = &vApp.stepdata[curStep].definitions;
        curIt = defMap->end();
        fillTable();
    }else if(change & (GUI::Change::definitions)){
        fillTable();
    }else if(change & (GUI::Change::atoms | GUI::Change::settings)){
        for(auto& def: *defMap){
            std::get<2>(def.second)->update(&std::get<0>(def.second),
                                            vApp.config.settings.atRadVdW.val,
                                            vApp.config.settings.atRadFac.val);
        }
    }
}

Vipster::Step::selection& DefineWidget::curSel()
{
    return std::get<0>(curIt->second);
}

Vipster::SelectionFilter& DefineWidget::curFilter()
{
    return std::get<1>(curIt->second);
}

std::shared_ptr<Vipster::GUI::SelData>& DefineWidget::curSelData()
{
    return std::get<2>(curIt->second);
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
    for(auto& [name, def]: *defMap){
        // if contained in extras, this group is displayed
        table.setItem(i, 0, new QTableWidgetItem{});
        table.item(i, 0)->setFlags(Qt::ItemIsSelectable|
                                   Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        table.item(i, 0)->setCheckState(Qt::CheckState(
                        master->curVP->hasExtraData(std::get<2>(def), false)*2));
        // name
        table.setItem(i, 1, new QTableWidgetItem(QString::fromStdString(name)));
        // filter-str
        table.setItem(i, 2, new QTableWidgetItem(QString::fromStdString(
            std::get<1>(def))));
        // color button
        auto* but = new QPushButton("Select");
        const auto& color = std::get<2>(def)->color;
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
        auto [it, _] = defMap->insert_or_assign(name,
            std::tuple{std::move(sel),
                filter,
                std::make_shared<GUI::SelData>()});
        curIt = it;
        curSelData()->update(&curSel(),
            vApp.config.settings.atRadVdW.val, vApp.config.settings.atRadFac.val);
        curSelData()->color = defaultColors[defMap->size()%5];
        master->curVP->addExtraData(curSelData(), false);
        triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
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
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: \"delete group\" triggered with invalid selection"};
    }
    master->curVP->delExtraData(curSelData(), false);
    defMap->erase(curIt);
    triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
}

void DefineWidget::on_fromSelButton_clicked()
{
    bool ok{false};
    auto tmp = QInputDialog::getText(this, "Copy selection to filtered group",
                                     "Enter name for new group:",
                                     QLineEdit::Normal, QString(), &ok
                                    ).toStdString();
    if(!ok) return;
    // convert selection to index filter
    const auto& idx = vApp.curSel->getAtoms().indices;
    std::string filter;
    std::stringstream ss{filter};
    ss << "[ ";
    for(const auto& p: idx){
        ss << p.first << " ";
    }
    ss << ']';
    // create new group
    auto [it, _] = defMap->insert_or_assign(tmp,
        std::tuple{*vApp.curSel,
                   filter,
                   std::make_shared<GUI::SelData>()});
    curIt = it;
    curSelData()->update(&curSel(),
        vApp.config.settings.atRadVdW.val, vApp.config.settings.atRadFac.val);
    curSelData()->color = defaultColors[defMap->size()%5];
    master->curVP->addExtraData(curSelData(), false);
    triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
}

void DefineWidget::toSelAction()
{
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: \"to selection\" triggered with invalid selection"};
    }
    *vApp.curSel = curSel();
    master->curVP->delExtraData(curSelData(), false);
    triggerUpdate(GUI::Change::selection | GUI::Change::extra);
}

void DefineWidget::updateAction()
{
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: \"update group\" triggered with invalid selection"};
    }
    curSel() = curStep->select(curFilter());
    curSelData()->update(&curSel(),
        vApp.config.settings.atRadVdW.val, vApp.config.settings.atRadFac.val);
    triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
}

void DefineWidget::on_defTable_cellChanged(int row, int column)
{
    if(column == 0){
        // checkbox consumes mouse event, need to select row manually
        ui->defTable->selectRow(row);
    }
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: selection is invalid."};
    }
    QTableWidgetItem *cell = ui->defTable->item(row, column);
    switch(column){
    case 0:
        // toggle visibility
        if(cell->checkState()){
            master->curVP->addExtraData(curSelData(), false);
        }else{
            master->curVP->delExtraData(curSelData(), false);
        }
        triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
        break;
    case 1:
        // change name
        if(defMap->find(cell->text().toStdString()) != defMap->end()){
            QMessageBox msg{this};
            msg.setText("Name \""+cell->text()+"\" is already in use.");
            msg.exec();
            QSignalBlocker block{ui->defTable};
            cell->setText(curIt->first.c_str());
            return;
        }else{
            auto node = defMap->extract(curIt);
            node.key() = cell->text().toStdString();
            defMap->insert(std::move(node));
        }
        break;
    case 2:
        // change filter
        try{
            auto filter = cell->text().toStdString();
            curSel() = curStep->select(filter);
            curSelData()->update(&curSel(),
                vApp.config.settings.atRadVdW.val, vApp.config.settings.atRadFac.val);
            triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
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
        curIt = defMap->end();
        for(auto* a: contextActions){
            a->setDisabled(true);
        }
    }else{
        const auto& curName = ui->defTable->item(sel[0]->row(), 1)->text();
        curIt = defMap->find(curName.toStdString());
        for(auto* a: contextActions){
            a->setEnabled(true);
        }
    }
}

void DefineWidget::colButton_clicked()
{
    auto& col = curSelData()->color;
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
    triggerUpdate(GUI::Change::definitions | GUI::Change::extra);
}

void DefineWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString{"About filters"}, FilterAbout);
}
