#include "pwparam.h"
#include "ui_pwparam.h"

using namespace Vipster;

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
    curParam = dynamic_cast<IO::PWParam*>(p);
    if(!curParam){
        throw Error("Invalid parameter set");
    }
    IO::PWParam::Namelist IO::PWParam::* namelists[] = {
        &IO::PWParam::control, &IO::PWParam::system,
        &IO::PWParam::electrons, &IO::PWParam::ions,
        &IO::PWParam::cell};
    for (int i=0; i<5; ++i) {
        auto* treeTop = ui->paramTree->topLevelItem(i);
        auto& nl = curParam->*namelists[i];
        // En-/Disable ions/cell namelists
        if(i>=3){
            if(nl.size()){
                treeTop->setFlags(treeTop->flags() | Qt::ItemIsEnabled);
            } else {
                treeTop->setFlags(treeTop->flags() & ~Qt::ItemIsEnabled);
            }
        }
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
}

void PWParam::addElement()
{
    auto& nl = curParam->*curNL;
    uint8_t idx{0};
    std::string key = "newKey";
    if(nl.find(key) != nl.end()){
        std::string key2 = key + "0";
        while(nl.find(key2) != nl.end()){
            key2 = key + std::to_string(++idx);
        }
        key = key2;
    }
    auto *child = new QTreeWidgetItem{{key.c_str(), nl[key].c_str()}};
    child->setFlags(child->flags() | Qt::ItemIsEditable);
    if(!curItem->parent()){
        curItem->addChild(child);
    }else{
        curItem->parent()->addChild(child);
    }
}

void PWParam::delElement()
{
    auto& nl = curParam->*curNL;
    nl.erase(nl.find(curKey));
    delete curItem;
}

void PWParam::on_paramTree_currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *)
{
    curItem = current;
    const std::map<QString, IO::PWParam::Namelist IO::PWParam::*> stringToNL =
    {
        {"control", &IO::PWParam::control},
        {"system", &IO::PWParam::system},
        {"electrons", &IO::PWParam::electrons},
        {"ions", &IO::PWParam::ions},
        {"cell", &IO::PWParam::cell},
    };
    if(!curItem->parent()){
        delAction->setDisabled(true);
        curNL = stringToNL.at(curItem->data(0, 0).toString().toLower());
    }else{
        delAction->setEnabled(true);
        curNL = stringToNL.at(curItem->parent()->data(0, 0).toString().toLower());
        curKey = curItem->data(0, 0).toString().toStdString();
    }
}

void PWParam::on_paramTree_itemChanged(QTreeWidgetItem *item, int column)
{
    auto& nl = curParam->*curNL;
    if(column == 0){
        auto newKey = item->data(0, 0).toString().toStdString();
        if(nl.find(newKey) != nl.end()){
            //TODO: alert user that he's trying to insert key twice?
            item->setData(0, 0, QString::fromStdString(curKey));
            return;
        }
        nl[newKey] = nl.at(curKey);
        nl.erase(nl.find(curKey));
        curKey = std::move(newKey);
    }else{
        nl.at(curKey) = item->data(1, 0).toString().toStdString();
    }
}
