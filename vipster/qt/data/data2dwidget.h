#ifndef DATA2DWIDGET_H
#define DATA2DWIDGET_H

#include <QWidget>
#include "../datawidget.h"
#include "../../common/meshdata.h"

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
    void on_sliceBut_toggled(bool checked);

private:
    Ui::Data2DWidget *ui;
    struct DatPlane{
        bool display;
        Vipster::GUI::MeshData gpu_data;
    };
    const Vipster::DataGrid2D_f* curData{nullptr};
    DatPlane* curPlane{nullptr};
    std::map<const Vipster::DataGrid2D_f*, DatPlane> planes;
};

#endif // DATA2DWIDGET_H
