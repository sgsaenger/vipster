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
    Vipster::BaseConfig* getConfig();
    Vipster::BaseParam* getParam();

private slots:
    void selFmt(int);

private:
    void enableParamWidget(bool);
    void enableConfWidget(bool);
    std::vector<Vipster::IOFmt> outFormats;
    Ui::SaveFmtDialog *ui;
};

#endif // SAVEFMTDIALOG_H
