#include "pwparam.h"
#include "ui_pwparam.h"

using namespace Vipster;

PWParam::PWParam(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PWParam)
{
    ui->setupUi(this);
}

PWParam::~PWParam()
{
    delete ui;
}

void PWParam::setParam(IO::PWParam* p)
{
    curParam = p;
    IO::PWNamelist IO::PWParam::* namelists[] = {
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
                treeTop->setFlags(treeTop->flags() & !Qt::ItemIsEnabled);
            }
        }
        // clear tree
        for (auto child: treeTop->takeChildren()) {
            delete child;
        }
        // add new entries
        for(auto& entry: nl){
            treeTop->addChild(new QTreeWidgetItem{{entry.first.c_str(), entry.second.c_str()}});
        }
    }
}
