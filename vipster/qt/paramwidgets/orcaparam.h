#ifndef ORCAPARAM_H
#define ORCAPARAM_H

#include <QWidget>
#include "../parambase.h"
#include "io/plugins/orca.h"

namespace Ui {
class ORCAParam;
}

class ORCAParam : public ParamBase
{
    Q_OBJECT

public:
    explicit ORCAParam(QWidget *parent = nullptr);
    ~ORCAParam();
    void setParam(Vipster::IO::Parameter *p) override;

private slots:
    void on_plainTextEdit_textChanged();

private:
    void saveText();
    Ui::ORCAParam *ui;
    Vipster::IO::Parameter* curParam{nullptr};
};

#endif // ORCAPARAM_H
