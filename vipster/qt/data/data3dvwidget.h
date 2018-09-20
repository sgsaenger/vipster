#ifndef DATA3DVWIDGET_H
#define DATA3DVWIDGET_H

#include <QWidget>
#include "../datawidget.h"

namespace Ui {
class Data3DVWidget;
}

class Data3DVWidget : public DataBase
{
    Q_OBJECT

public:
    explicit Data3DVWidget(QWidget *parent = nullptr);
    ~Data3DVWidget() override;
    void setData(const Vipster::BaseData *data) override;

private:
    const Vipster::DataGrid3D_v* curData{nullptr};
    Ui::Data3DVWidget *ui;
};

#endif // DATA3DVWIDGET_H
