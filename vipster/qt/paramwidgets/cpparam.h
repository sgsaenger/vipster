#ifndef CPPARAM_H
#define CPPARAM_H

#include <QWidget>
#include "../parambase.h"
#include "io/plugins/cpmdinput.h"

namespace Ui {
class CPParam;
}

class CPParam : public ParamBase
{
    Q_OBJECT

public:
    explicit CPParam(QWidget *parent = nullptr);
    ~CPParam() override;
    void setParam(Vipster::IO::Parameter *p) override;

private slots:
    void on_comboBox_currentIndexChanged(const QString &arg);
    void on_plainTextEdit_textChanged();
    void on_prefixEdit_editingFinished();
    void on_suffixEdit_editingFinished();
    void on_nlEdit_editingFinished();

private:
    void fillText();
    void saveText();
    Ui::CPParam *ui;
    Vipster::IO::Parameter* curParam{nullptr};
    std::vector<std::string>* curSection{nullptr};
};

#endif // CPPARAM_H
