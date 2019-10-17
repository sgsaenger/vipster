#ifndef MILLERWIDGET_H
#define MILLERWIDGET_H

#include <QWidget>
#include "../basewidget.h"
#include "../../common/meshdata.h"
#include "molecule.h"

namespace Ui {
class MillerWidget;
}

class MillerWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MillerWidget(QWidget *parent = nullptr);
    ~MillerWidget() override;
    void updateWidget(Vipster::GUI::change_t) override;

private slots:
    void updateIndex(int idx);
    void updateOffset(double off);
    void on_pushButton_toggled(bool checked);

private:
    Ui::MillerWidget *ui;
    struct MillerPlane: Vipster::GUI::MeshData
    {
        MillerPlane(const Vipster::GUI::GlobalData& glob, std::vector<Face>&& faces,
                    Vipster::Vec offset, Vipster::Mat cell, Texture texture,
                    const std::array<int8_t, 3> &hkl);
        std::array<int8_t, 3> hkl;
    };
    Vipster::Step* curStep{nullptr};
    MillerPlane* curPlane{nullptr};
};

#endif // MILLERWIDGET_H
