#ifndef CONFIGWIDGET_H
#define CONFIGWIDGET_H

#include <QWidget>
#include "basewidget.h"
#include "io.h"

namespace Ui {
class PresetWidget;
}

class PresetBase: public QWidget
{
    Q_OBJECT

public:
    explicit PresetBase(QWidget *parent = nullptr);
    virtual ~PresetBase() = default;
    virtual void setPreset(Vipster::IO::BasePreset *c)=0;
};

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
    Vipster::IOFmt curFmt;
    Vipster::IO::BasePreset* curPreset{nullptr};

private slots:
    void on_presetSel_currentIndexChanged(int index);

    void on_helpButton_clicked();

private:
    Ui::PresetWidget *ui;
};

#endif // CONFIGWIDGET_H
