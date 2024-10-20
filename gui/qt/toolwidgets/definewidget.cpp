#include "definewidget.h"
#include "ui_definewidget.h"
#include "mainwindow.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QColorDialog>
#include <QInputDialog>

using namespace Vipster;

DefineWidget::DefineWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DefineWidget)
{
    ui->setupUi(this);
    contextActions.push_back(new QAction{"Update group", ui->defTable});
    contextActions.back()->setDisabled(true);
    connect(contextActions.back(), &QAction::triggered, this, &DefineWidget::updateDefinition);
    contextActions.push_back(new QAction{"Delete group", ui->defTable});
    contextActions.back()->setDisabled(true);
    connect(contextActions.back(), &QAction::triggered, this, &DefineWidget::deleteDefinition);
    contextActions.push_back(new QAction{"Set as selection", ui->defTable});
    contextActions.back()->setDisabled(true);
    connect(contextActions.back(), &QAction::triggered, this, &DefineWidget::copyDefToSelection);
    ui->defTable->addActions(contextActions);

    // load definition map for current step
    connect(&vApp, &Application::activeStepChanged,
            this, [&](const Step &step) {
        defMap = &vApp.getState(step).definitions;
        curIt = defMap->end();
        fillTable();
    });

    // update definitions if required
    connect(&vApp, &Application::stepChanged,
            this, [&](const Step &step) {
        if (&step != &vApp.curStep()) return;
        const auto &settings = vApp.config().settings;
        for (auto &[name, def]: *defMap) {
            auto &[selection, filter, selData] = def;
            selData->update(&selection,
                            settings.atRadVdW.val,
                            settings.atRadFac.val);
        }
    });

    // update visualization settings
    // TODO: this only updates the current step's define's visualization.
    // TODO: maybe the renderData should not be persistent but created on the fly? less state!
    connect(&vApp, &Application::configChanged,
            this, [&](const ConfigState &c){
        const auto &settings = c.settings;
        for (auto &[name, def]: *defMap) {
            auto &[selection, filter, selData] = def;
            selData->update(&selection,
                            settings.atRadVdW.val,
                            settings.atRadFac.val);
        }
    });

    // TODO: scriptwidget may also create definitions, how to interact?

    connect(ui->newButton, &QPushButton::clicked, this, &DefineWidget::createDefinition);
    connect(ui->fromSelButton, &QPushButton::clicked, this, &DefineWidget::copySelToDefinition);
    connect(ui->helpButton, &QPushButton::clicked,
            this, [&](){ QMessageBox::information(this, QString{"About filters"}, FilterAbout); });

    connect(ui->defTable, &QTableWidget::cellChanged, this, &DefineWidget::tableCellChanged);
    connect(ui->defTable, &QTableWidget::itemSelectionChanged, this, &DefineWidget::tableSelectionChanged);
}

DefineWidget::~DefineWidget()
{
    delete ui;
    for(auto* i: contextActions){
        delete i;
    }
}

Vipster::Step::const_selection& DefineWidget::curSel()
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
    for(const auto &[name, def]: *defMap){
        const auto &[sel, filter, selData] = def;
        // if contained in extras, this group is displayed
        table.setItem(i, 0, new QTableWidgetItem{});
        table.item(i, 0)->setFlags(Qt::ItemIsSelectable|
                                   Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
        table.item(i, 0)->setCheckState(Qt::CheckState(
            static_cast<MainWindow*>(this->parent()->parent())->curVP->hasExtraData(selData, false)*2));
        // name
        table.setItem(i, 1, new QTableWidgetItem(QString::fromStdString(name)));
        // filter-str
        table.setItem(i, 2, new QTableWidgetItem(QString::fromStdString(filter)));
        // color button
        auto* but = new QPushButton("Select");
        const auto& color = selData->color;
        but->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                           .arg(color[0]).arg(color[1]).arg(color[2]));
        connect(but, &QPushButton::clicked, this, &DefineWidget::changeColor);
        table.setCellWidget(i, 3, but);
        i++;
    }
}

void DefineWidget::createDefinition()
{
    bool ok{false};
    const auto filter = QInputDialog::getText(this, "Create new filtered group",
                                              "Enter filter for new group:",
                                              QLineEdit::Normal, QString(), &ok
                                             ).toStdString();
    if(!ok) return;
    try {
        auto sel = vApp.curStep().select(filter);

        const auto name = QInputDialog::getText(this, "Create new filtered group",
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
            vApp.config().settings.atRadVdW.val, vApp.config().settings.atRadFac.val);
        curSelData()->color = defaultColors[defMap->size()%5];

        // TODO: find better access method. this widget is reparented into a dock widget -> implicit two layers of indirection
        static_cast<MainWindow*>(this->parent()->parent())->curVP->addExtraData(curSelData(), false);
        fillTable();
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

void DefineWidget::deleteDefinition()
{
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: \"delete group\" triggered with invalid selection"};
    }
    static_cast<MainWindow*>(this->parent()->parent())->curVP->delExtraData(curSelData(), false);
    defMap->erase(curIt);
    fillTable();
}

void DefineWidget::copySelToDefinition()
{
    bool ok{false};
    auto tmp = QInputDialog::getText(this, "Copy selection to filtered group",
                                     "Enter name for new group:",
                                     QLineEdit::Normal, QString(), &ok
                                    ).toStdString();
    if(!ok) return;

    // convert selection to index filter
    const auto& idx = vApp.curSel().getAtoms().indices;
    std::stringstream ss{};
    ss << "index [ ";
    for(const auto& p: idx){
        ss << p.first << " ";
    }
    ss << ']';

    // create new group
    Step::const_selection constSel = vApp.curSel();
    auto [it, _] = defMap->insert_or_assign(tmp,
        std::tuple{constSel,
                   ss.str(),
                   std::make_shared<GUI::SelData>()});
    curIt = it;
    curSelData()->update(&curSel(),
        vApp.config().settings.atRadVdW.val, vApp.config().settings.atRadFac.val);
    curSelData()->color = defaultColors[defMap->size()%5];

    static_cast<MainWindow*>(this->parent()->parent())->curVP->addExtraData(curSelData(), false);
    fillTable();
}

void DefineWidget::copyDefToSelection()
{
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: \"to selection\" triggered with invalid selection"};
    }
    // copy selection to definition
    vApp.updateSelection(curFilter());
    // hide definition
    static_cast<MainWindow*>(this->parent()->parent())->curVP->delExtraData(curSelData(), false);
}

void DefineWidget::updateDefinition()
{
    if(curIt == defMap->end()){
        throw Error{"DefineWidget: \"update group\" triggered with invalid selection"};
    }
    curSel() = vApp.curStep().select(curFilter());
    curSelData()->update(&curSel(),
        vApp.config().settings.atRadVdW.val, vApp.config().settings.atRadFac.val);
    static_cast<MainWindow*>(this->parent()->parent())->curVP->updateState();
}

void DefineWidget::tableCellChanged(int row, int column)
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
            static_cast<MainWindow*>(this->parent()->parent())->curVP->addExtraData(curSelData(), false);
        }else{
            static_cast<MainWindow*>(this->parent()->parent())->curVP->delExtraData(curSelData(), false);
        }
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
            curFilter() = filter;
            curSel() = vApp.curStep().select(filter);
            curSelData()->update(&curSel(),
                vApp.config().settings.atRadVdW.val, vApp.config().settings.atRadFac.val);
            static_cast<MainWindow*>(this->parent()->parent())->curVP->updateState();
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

void DefineWidget::tableSelectionChanged()
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

void DefineWidget::changeColor()
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
    static_cast<MainWindow*>(this->parent()->parent())->curVP->updateState();
}
