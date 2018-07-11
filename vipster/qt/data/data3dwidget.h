#ifndef DATA3DWIDGET_H
#define DATA3DWIDGET_H

#include <QWidget>
#include "../datawidget.h"

namespace Ui {
class Data3DWidget;
}

class Data3DWidget : public DataBase
{
    Q_OBJECT

public:
    explicit Data3DWidget(QWidget *parent = nullptr);
    ~Data3DWidget() override;
    void setData(const Vipster::BaseData *data) override;

private slots:
    void on_sliceBut_clicked();
    void on_sliceDir_currentIndexChanged(int index);
    void on_sliceVal_valueChanged(int value);
    void on_surfToggle_stateChanged(int arg1);
    void on_surfVal_valueChanged(double arg1);
    void on_surfBut_clicked();

private:
    const Vipster::DataGrid3D_f* curData{nullptr};
    Ui::Data3DWidget *ui;
};

#endif // DATA3DWIDGET_H
