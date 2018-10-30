#ifndef DATA3DWIDGET_H
#define DATA3DWIDGET_H

#include <QWidget>
#include <QDoubleValidator>
#include "../datawidget.h"
#include "../../common/meshdata.h"

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
    void updateWidget(Vipster::guiChange_t change) override;

private slots:
    void on_sliceDir_currentIndexChanged(int index);
    void on_sliceVal_valueChanged(int value);
    void on_surfToggle_stateChanged(int arg1);
    void on_surfSlider_valueChanged(int value);
    void on_surfVal_editingFinished();
    void on_sliceBut_toggled(bool checked);
    void on_surfBut_toggled(bool checked);

private:
    Ui::Data3DWidget *ui;
    const Vipster::DataGrid3D_f* curData{nullptr};
    QDoubleValidator validator;

    struct DatSlice{
        bool display;
        size_t dir, pos;
        Vipster::GUI::MeshData gpu_data;
    };
    DatSlice* curSlice{nullptr};
    std::map<const Vipster::DataGrid3D_f*, DatSlice> slices;

    struct IsoSurf{
        bool display, plusmin;
        float isoval;
        Vipster::GUI::MeshData gpu_data;
    };
    IsoSurf* curSurf{nullptr};
    std::map<const Vipster::DataGrid3D_f*, IsoSurf> surfaces;
    std::vector<Vipster::GUI::MeshData::Face> mkSurf(float isoval, bool pm);
};

#endif // DATA3DWIDGET_H
