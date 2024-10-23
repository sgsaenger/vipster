#ifndef DATA2DWIDGET_H
#define DATA2DWIDGET_H

#include <memory>
#include <QWidget>
#include "../mainwidgets/datawidget.h"
#include "meshdata.h"

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
    const Vipster::DataGrid2D_f* curData{nullptr};
    std::map<const Vipster::DataGrid2D_f*,
             std::shared_ptr<Vipster::GUI::MeshData>> planes;
    std::shared_ptr<Vipster::GUI::MeshData> curPlane{nullptr};
};

#endif // DATA2DWIDGET_H
