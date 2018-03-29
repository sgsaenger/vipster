#ifndef SAVEFMTDIALOG_H
#define SAVEFMTDIALOG_H

#include <QDialog>
#include <vector>
#include "iowrapper.h"

namespace Ui {
class SaveFmtDialog;
}

class SaveFmtDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SaveFmtDialog(QWidget *parent = nullptr);
    ~SaveFmtDialog();
    Vipster::IOFmt fmt{};
    Vipster::BaseParam *param{nullptr};
    Vipster::BaseConfig *config{nullptr};

private slots:
    void selFmt(int);
    void on_paramSel_currentIndexChanged(int index);
    void on_configSel_currentIndexChanged(int index);

private:
    void enableParamWidget(bool);
    void enableConfWidget(bool);
    std::vector<Vipster::IOFmt> outFormats;
    std::vector<std::unique_ptr<Vipster::BaseParam>> ownParams;
    std::vector<std::unique_ptr<Vipster::BaseConfig>> ownConfigs;
    Ui::SaveFmtDialog *ui;
    static constexpr int paramlist[] = {0, 1, 0, 0, 0};
    static constexpr int conflist[] =  {0, 0, 0, 0, 0};
};

#endif // SAVEFMTDIALOG_H
