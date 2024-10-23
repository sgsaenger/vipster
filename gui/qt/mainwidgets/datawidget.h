#ifndef DATAWIDGET_H
#define DATAWIDGET_H

#include <QWidget>
#include "vipster/data.h"

namespace Ui {
class DataWidget;
}

class DataBase: public QWidget
{
    Q_OBJECT

public:
    explicit DataBase(QWidget *parent = nullptr);
    virtual void setData(const Vipster::BaseData *d)=0;
};

class DataWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DataWidget(QWidget *parent = nullptr);
    ~DataWidget() override;

private slots:
    void selectData(int index);

private:
    Ui::DataWidget *ui;
    const Vipster::BaseData* curData{nullptr};
};

#endif // DATAWIDGET_H
