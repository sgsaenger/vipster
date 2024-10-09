#ifndef CONFIGWIDGET_H
#define CONFIGWIDGET_H

#include <QWidget>
#include <QHBoxLayout>
#include "vipster/fileio.h"

namespace Ui {
class PresetWidget;
}

class PresetWidget : public QWidget
{
    Q_OBJECT

public:
    explicit PresetWidget(QWidget *parent = nullptr);
    ~PresetWidget() override;
    std::vector<std::pair<std::string, Vipster::Preset>> presets;
    void registerPreset(const std::string& name,
                        const Vipster::Preset& data);
    void clearPresets();
    Vipster::Preset* curPreset{nullptr};

private slots:
    void on_presetSel_currentIndexChanged(int index);

    void on_helpButton_clicked();

private:
    Ui::PresetWidget *ui;
};

#endif // CONFIGWIDGET_H
