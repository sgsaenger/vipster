#ifndef PARAMWIDGET_H
#define PARAMWIDGET_H

#include "../basewidget.h"
#include "../paramwidgets.h"
#include "io.h"

namespace Ui {
class ParamWidget;
}

class ParamWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit ParamWidget(QWidget *parent = nullptr);
    ~ParamWidget() override;
    std::vector<std::pair< Vipster::IOFmt, std::unique_ptr<Vipster::IO::BaseParam>>> params;
    void registerParam(std::unique_ptr<Vipster::IO::BaseParam>&& data);
    void clearParams();
    Vipster::IOFmt curFmt;
    Vipster::IO::BaseParam *curParam{nullptr};

private slots:
    void on_paramSel_currentIndexChanged(int index);

    void on_pushButton_clicked();

private:
    Ui::ParamWidget *ui;
    std::map<Vipster::IOFmt, ParamBase*> formats;
};

#endif // PARAMWIDGET_H
