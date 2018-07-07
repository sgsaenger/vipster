#ifndef PWIPARAM_H
#define PWIPARAM_H

#include <QWidget>
#include <QTreeWidgetItem>
#include <QAction>
#include "io/pwinput/param.h"
#include "../paramwidget.h"

namespace Ui {
class PWParam;
}

class PWParam : public ParamBase
{
    Q_OBJECT

public:
    explicit PWParam(QWidget *parent = nullptr);
    ~PWParam() override;
    void setParam(Vipster::IO::BaseParam *p) override;

private slots:
    void addElement();
    void delElement();
    void on_paramTree_currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *previous);
    void on_paramTree_itemChanged(QTreeWidgetItem *item, int column);

private:
    Ui::PWParam *ui;
    QAction *addAction;
    QAction *delAction;
    QTreeWidgetItem *curItem{nullptr};
    std::string curKey{};
    Vipster::IO::PWParam::Namelist Vipster::IO::PWParam::* curNL{nullptr};
    Vipster::IO::PWParam *curParam{nullptr};
};

#endif // PWIPARAM_H
