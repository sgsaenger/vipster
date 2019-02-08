#ifndef ORCAPARAM_H
#define ORCAPARAM_H

#include <QWidget>
#include "../paramwidget.h"
#include "io/orca/param.h"

namespace Ui {
class ORCAParam;
}

class ORCAParam : public ParamBase
{
    Q_OBJECT

public:
    explicit ORCAParam(QWidget *parent = nullptr);
    ~ORCAParam();
    void setParam(Vipster::IO::BaseParam *p) override;

private slots:
    void on_plainTextEdit_textChanged();

private:
    void saveText();
    Ui::ORCAParam *ui;
    Vipster::IO::OrcaParam* curParam{nullptr};
};

#endif // ORCAPARAM_H
