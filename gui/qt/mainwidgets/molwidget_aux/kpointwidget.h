#ifndef KPOINTWIDGET_H
#define KPOINTWIDGET_H

#include <QWidget>
#include "vipster/molecule.h"

namespace Ui {
    class KPointWidget;
}

class KPointWidget : public QWidget
{
    Q_OBJECT

public:
    explicit KPointWidget(QWidget *parent = nullptr);
    ~KPointWidget();

private:
    void setActiveMol(const Vipster::Molecule &mol);
    void updateMol(const Vipster::Molecule &mol);

    Ui::KPointWidget *ui;
};

#endif // KPOINTWIDGET_H
