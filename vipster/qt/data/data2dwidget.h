#ifndef DATA2DWIDGET_H
#define DATA2DWIDGET_H

#include <QWidget>
#include "../datawidget.h"

namespace Ui {
class Data2DWidget;
}

class Data2DWidget : public DataBase
{
    Q_OBJECT

public:
    explicit Data2DWidget(QWidget *parent = nullptr);
    ~Data2DWidget() override;
    void setData(const Vipster::BaseData *d) override;

private slots:
    void on_sliceBut_clicked();
    void on_sliceDir_currentIndexChanged(int index);
    void on_sliceVal_valueChanged(int value);

private:
    const Vipster::DataGrid2D_f* curData{nullptr};
    Ui::Data2DWidget *ui;
};

#endif // DATA2DWIDGET_H
