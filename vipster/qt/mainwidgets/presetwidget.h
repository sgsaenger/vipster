#ifndef CONFIGWIDGET_H
#define CONFIGWIDGET_H

#include <QWidget>
#include "../basewidget.h"
#include "../presetwidgets.h"
#include "io.h"

namespace Ui {
class PresetWidget;
}

class PresetWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PresetWidget(QWidget *parent = nullptr);
    ~PresetWidget() override;
    std::vector<std::pair<std::string, std::unique_ptr<Vipster::IO::BasePreset>>> presets;
    void registerPreset(const std::string& name,
                        std::unique_ptr<Vipster::IO::BasePreset>&& data);
    void clearPresets();
    const Vipster::IO::Plugin* curFmt{nullptr};
    Vipster::IO::BasePreset* curPreset{nullptr};

private slots:
    void on_presetSel_currentIndexChanged(int index);

    void on_helpButton_clicked();

private:
    Ui::PresetWidget *ui;
    std::map<const Vipster::IO::Plugin*, PresetBase*> formats;
};

#endif // CONFIGWIDGET_H
