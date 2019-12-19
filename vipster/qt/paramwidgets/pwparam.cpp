#include "pwparam.h"
#include "ui_pwparam.h"
#include <QMessageBox>

using namespace Vipster;

using NameList = std::map<std::string, std::string>;
static const std::string namelists[]{
    "&CONTROL",
    "&SYSTEM",
    "&ELECTRONS",
    "&IONS",
    "&CELL"
};

PWParam::PWParam(QWidget *parent) :
    ParamBase(parent),
    ui(new Ui::PWParam)
{
    ui->setupUi(this);
    addAction = new QAction{"Add Element", ui->paramTree};
    ui->paramTree->addAction(addAction);
    connect(addAction, &QAction::triggered, this, &PWParam::addElement);
    delAction = new QAction{"Delete Element", ui->paramTree};
    ui->paramTree->addAction(delAction);
    connect(delAction, &QAction::triggered, this, &PWParam::delElement);
}

PWParam::~PWParam()
{
    delete ui;
    delete addAction;
    delete delAction;
}

void PWParam::setParam(IO::BaseParam *p)
{
    auto treeBlocker = QSignalBlocker{ui->paramTree};
    curParam = p;
    if(curParam->getFmt() != &IO::PWInput){
        throw Error("Invalid parameter set");
    }
    for (int i=0; i<5; ++i) {
        auto* treeTop = ui->paramTree->topLevelItem(i);
        auto& nl = std::get<NameList>(curParam->at(namelists[i]));
        // clear tree
        for (auto child: treeTop->takeChildren()) {
            delete child;
        }
        // add new entries
        for(auto& entry: nl){
            auto *child = new QTreeWidgetItem{{entry.first.c_str(), entry.second.c_str()}};
            child->setFlags(child->flags() | Qt::ItemIsEditable);
            treeTop->addChild(child);
        }
    }
    ui->prefixEdit->setText(std::get<std::string>(curParam->at("PPPrefix")).c_str());
    ui->suffixEdit->setText(std::get<std::string>(curParam->at("PPSuffix")).c_str());
}

void PWParam::addElement()
{
    std::string key = "newKey";
    if(curNL->find(key) != curNL->end()){
        uint8_t idx{0};
        std::string key2 = key + "0";
        while(curNL->find(key2) != curNL->end()){
            key2 = key + std::to_string(++idx);
        }
        key = key2;
    }
    auto *child = new QTreeWidgetItem{{key.c_str(), (*curNL)[key].c_str()}};
    child->setFlags(child->flags() | Qt::ItemIsEditable);
    if(!curItem->parent()){
        curItem->addChild(child);
    }else{
        curItem->parent()->addChild(child);
    }
}

void PWParam::delElement()
{
    curNL->erase(curNL->find(curKey));
    delete curItem;
}

void PWParam::on_paramTree_currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *)
{
    curItem = current;
    if(!curItem->parent()){
        delAction->setDisabled(true);
        curNL = &std::get<NameList>(curParam->at(curItem
                     ->data(0,0).toString().toStdString()));
    }else{
        delAction->setEnabled(true);
        curNL = &std::get<NameList>(curParam->at(curItem->parent()
                     ->data(0,0).toString().toStdString()));
        curKey = curItem->data(0, 0).toString().toStdString();
    }
}

void PWParam::on_paramTree_itemChanged(QTreeWidgetItem *item, int column)
{
    if(column == 0){
        auto newKey = item->data(0, 0).toString().toStdString();
        if(curNL->find(newKey) != curNL->end()){
            QMessageBox::warning(this, "Key already present",
                                 QString{"Cannot rename key \""}+curKey.c_str()+"\" to \""+
                                 newKey.c_str()+"\" because the latter is already present.");
            QSignalBlocker block{ui->paramTree};
            item->setData(0, 0, QString::fromStdString(curKey));
            return;
        }
        (*curNL)[newKey] = curNL->at(curKey);
        curNL->erase(curNL->find(curKey));
        curKey = std::move(newKey);
    }else{
        curNL->at(curKey) = item->data(1, 0).toString().toStdString();
    }
}

void PWParam::on_prefixEdit_editingFinished()
{
    std::get<std::string>(curParam->at("PPPrefix")) =
            ui->prefixEdit->text().toStdString();
}

void PWParam::on_suffixEdit_editingFinished()
{
    std::get<std::string>(curParam->at("PPSuffix")) =
            ui->suffixEdit->text().toStdString();
}
