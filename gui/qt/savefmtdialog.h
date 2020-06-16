#ifndef SAVEFMTDIALOG_H
#define SAVEFMTDIALOG_H

#include <QDialog>
#include <vector>
#include "vipster/fileio.h"

namespace Ui {
class SaveFmtDialog;
}

class SaveFmtDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SaveFmtDialog(const Vipster::PluginList& plugins, QWidget *parent = nullptr);
    ~SaveFmtDialog();
    const Vipster::Plugin* plugin{};
    std::optional<Vipster::Preset> getPreset();
    std::optional<Vipster::Parameter> getParam();

private slots:
    void selFmt(int);

private:
    void enableParamWidget(bool);
    void enablePresetWidget(bool);
    std::vector<const Vipster::Plugin*> outFormats;
    Ui::SaveFmtDialog *ui;
};

#endif // SAVEFMTDIALOG_H
