#ifndef DATA3DWIDGET_H
#define DATA3DWIDGET_H

#include <memory>
#include <QWidget>
#include <QDoubleValidator>
#include "../mainwidgets/datawidget.h"
#include "meshdata.h"

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
    void updateWidget(Vipster::GUI::change_t change) override;

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

    struct DatSlice: public Vipster::GUI::MeshData{
        size_t dir, pos;
        DatSlice(std::vector<Face>&& faces,
                 Vipster::Vec offset, Vipster::Mat cell,
                 Texture texture, size_t dir, size_t pos);
    };
    std::shared_ptr<DatSlice> curSlice{nullptr};
    std::map<const Vipster::DataGrid3D_f*, std::shared_ptr<DatSlice>> slices;

    struct IsoSurf: public Vipster::GUI::MeshData{
        bool plusmin;
        double isoval;
        IsoSurf(std::vector<Face>&& faces,
                Vipster::Vec offset, Vipster::Mat cell,
                Texture texture, bool plusmin, double isoval);
    };
    std::shared_ptr<IsoSurf> curSurf{nullptr};
    std::map<const Vipster::DataGrid3D_f*, std::shared_ptr<IsoSurf>> surfaces;
    std::vector<Vipster::GUI::MeshData::Face> mkSurf(double isoval, bool pm);
};

#endif // DATA3DWIDGET_H
