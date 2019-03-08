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

void DefineWidget::updateWidget(Vipster::guiChange_t change)
{
    if((change & guiStepChanged) == guiStepChanged){
        for(auto& pair: dataMap[curStep]){
            if(pair.second.display) master->delExtraData(&pair.second.gpu_data);
        }
        curStep = master->curStep;
        defMap = &master->stepdata[curStep].def;
        fillTable();
        for(auto& pair: dataMap[curStep]){
            if(pair.second.display) master->addExtraData(&pair.second.gpu_data);
        }
    }else if(change & GuiChange::definitions){
        fillTable();
    }else if(change & GuiChange::atoms){
        auto& curMap = dataMap[curStep];
        for(auto& name: curNames){
            curMap.at(name).gpu_data.update(&defMap->at(name));
        }
    }
}

void DefineWidget::fillTable()
{
    ui->defTable->clearSelection();
    QSignalBlocker blockTable{ui->defTable};
    auto& curMap = dataMap[curStep];
    // make sure data is up to date and synchronized
    for(auto& def: *defMap){
        const auto& name = def.first;
        def.second.evaluateCache();
        auto pos = curMap.find(name);
        if(pos == curMap.end()){
            auto& curCol = colors[curMap.size()%5];
            auto tmp = curMap.emplace(name, GroupData{true, curCol,
                                           GUI::SelData{master->getGLGlobals(),
                                                        curCol,
                                                        &def.second}});
            master->addExtraData(&tmp.first->second.gpu_data);
        }
    }
    for(auto it = curMap.begin(); it != curMap.end();){
        auto pos = defMap->find(it->first);
        if(pos == defMap->end()){
            if(it->second.display){
                master->delExtraData(&it->second.gpu_data);
            }
            it = curMap.erase(it);
        }else{
            it->second.gpu_data.update(&pos->second);
            ++it;
        }
    }
    // setup table
    int i{0};
    auto& table = *ui->defTable;
    table.clearContents();
    curNames.clear();
    ui->defTable->setRowCount(defMap->size());
    for(auto& pair: curMap){
        curNames.push_back(pair.first);
        auto& dat = pair.second;
        table.setItem(i, 0, new QTableWidgetItem(QString::fromStdString(pair.first)));
        table.item(i, 0)->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEditable|
                                   Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        table.item(i, 0)->setCheckState(Qt::CheckState(static_cast<int>(dat.display)*2));
        table.setItem(i, 1, new QTableWidgetItem(QString::fromStdString(
                                                     defMap->at(pair.first).getFilter().toStr())));
        auto* but = new QPushButton("Select");
        but->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                           .arg(dat.color[0]).arg(dat.color[1]).arg(dat.color[2]));
        connect(but, &QPushButton::clicked, this, &DefineWidget::colButton_clicked);
        table.setCellWidget(i, 2, but);
        i++;
    }
}

void DefineWidget::on_newButton_clicked()
{
    auto filter = QInputDialog::getText(this, "Create new filtered group",
                                        "Enter filter for new group:").toStdString();
    try {
        auto tmp = curStep->select(filter);
        auto name = QInputDialog::getText(this, "Create new filtered group",
                                          "Enter name for new group:").toStdString();
        defMap->insert_or_assign(name, std::move(tmp));
        triggerUpdate(GuiChange::definitions);
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
    if(curSel < 0){
        throw Error{"DefineWidget: \"delete group\" triggered with invalid selection"};
    }
    defMap->erase(curNames.at(curSel));
    triggerUpdate(GuiChange::definitions);
}

void DefineWidget::on_fromSelButton_clicked()
{
    auto tmp = QInputDialog::getText(this, "Copy selection to filtered group",
                                     "Enter name for new group:").toStdString();
    defMap->insert_or_assign(tmp, *master->curSel);
    triggerUpdate(GuiChange::definitions);
}

void DefineWidget::toSelAction()
{
    if(curSel < 0){
        throw Error{"DefineWidget: \"to selection\" triggered with invalid selection"};
    }
    *master->curSel = defMap->at(curNames.at(curSel));
    triggerUpdate(GuiChange::selection);
}

void DefineWidget::updateAction()
{
    if(curSel < 0){
        throw Error{"DefineWidget: \"update group\" triggered with invalid selection"};
    }
    auto& step = defMap->at(curNames[curSel]);
    step.setFilter(step.getFilter());
    triggerUpdate(GuiChange::definitions);
}

void DefineWidget::on_defTable_cellChanged(int row, int column)
{
    QTableWidgetItem *cell = ui->defTable->item(row, column);
    auto& name = curNames[static_cast<size_t>(row)];
    auto& curMap = dataMap[curStep];
    auto& curData = curMap.at(name);
    if(column == 0){
        bool checkState = cell->checkState() != Qt::CheckState::Unchecked;
        if(checkState != curData.display){
            curData.display = checkState;
            if(checkState){
                master->addExtraData(&curData.gpu_data);
            }else{
                master->delExtraData(&curData.gpu_data);
            }
        }else{
            auto newName = cell->text();
            if(defMap->find(newName.toStdString()) != defMap->end()){
                QSignalBlocker block{ui->defTable};
                cell->setText(name.c_str());
                QMessageBox msg{this};
                msg.setText("Name \""+newName+"\" is already in use.");
                msg.exec();
                return;
            }
            auto it1 = defMap->find(name);
            auto it2 = curMap.find(name);
            name = newName.toStdString();
            master->delExtraData(&it2->second.gpu_data);
            defMap->emplace(name, std::move(it1->second));
            defMap->erase(it1);
            auto res = curMap.emplace(name, std::move(it2->second));
            curMap.erase(it2);
            master->addExtraData(&res.first->second.gpu_data);
            // TODO: reenable when libc++-8 is available
//            auto nh1 = defMap->extract(name);
//            auto nh2 = curMap.extract(name);
//            name = cell->text().toStdString();
//            nh1.key() = name;
//            nh2.key() = name;
//            defMap->insert(std::move(nh1));
//            curMap.insert(std::move(nh2));
        }
    }else{
        try{
            auto filter = cell->text().toStdString();
            defMap->at(name).setFilter(filter);
            triggerUpdate(GuiChange::definitions);
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
}

void DefineWidget::on_defTable_itemSelectionChanged()
{
    auto sel = ui->defTable->selectedItems();
    if(sel.empty()){
        curSel = -1;
        for(auto* a: contextActions){
            a->setDisabled(true);
        }
    }else{
        curSel = sel[0]->row();
        for(auto* a: contextActions){
            a->setEnabled(true);
        }
    }
}

void DefineWidget::colButton_clicked()
{
    auto& curData = dataMap.at(curStep).at(curNames.at(curSel));
    auto& col = curData.color;
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
    curData.gpu_data.update(col);
    static_cast<QPushButton*>(sender())->setStyleSheet(
        QString("background-color: %1").arg(newCol.name()));
    triggerUpdate(GuiChange::definitions);
}

void DefineWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString{"About filters"}, FilterAbout);
}
