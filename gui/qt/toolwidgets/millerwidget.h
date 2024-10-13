#ifndef MILLERWIDGET_H
#define MILLERWIDGET_H

#include <QWidget>
#include "meshdata.h"
#include "mainwindow.h"
#include "vipster/molecule.h"

namespace Ui {
class MillerWidget;
}

class MillerWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MillerWidget(QWidget *parent);
    ~MillerWidget() override;

private slots:
    void updateIndex(int idx);
    void updateOffset(double off);

private:
    Ui::MillerWidget *ui;
    struct MillerPlane: Vipster::GUI::MeshData
    {
        MillerPlane(std::vector<Face>&& faces, Vipster::Vec offset,
                    Vipster::Mat cell, Texture texture,
                    const std::array<int8_t, 3> &hkl);
        std::array<int8_t, 3> hkl;
    };
    // TODO: use weak_ptr
    std::map<const Vipster::Step*, std::shared_ptr<MillerPlane>> planes;
    std::shared_ptr<MillerPlane> curPlane{nullptr};
};

#endif // MILLERWIDGET_H
