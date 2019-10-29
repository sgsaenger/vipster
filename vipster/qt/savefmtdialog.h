#ifndef SAVEFMTDIALOG_H
#define SAVEFMTDIALOG_H

#include <QDialog>
#include <vector>
#include "io.h"

namespace Ui {
class SaveFmtDialog;
}

class SaveFmtDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SaveFmtDialog(const Vipster::IO::Plugins& plugins, QWidget *parent = nullptr);
    ~SaveFmtDialog();
    const Vipster::IO::Plugin* plugin{};
    Vipster::IO::BaseConfig* getConfig();
    Vipster::IO::BaseParam* getParam();

private slots:
    void selFmt(int);

private:
    void enableParamWidget(bool);
    void enableConfWidget(bool);
    std::vector<const Vipster::IO::Plugin*> outFormats;
    Ui::SaveFmtDialog *ui;
};

#endif // SAVEFMTDIALOG_H
