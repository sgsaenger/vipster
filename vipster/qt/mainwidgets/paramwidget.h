#ifndef PARAMWIDGET_H
#define PARAMWIDGET_H

#include "../basewidget.h"
#include "../paramwidgets.h"
#include "fileio.h"

namespace Ui {
class ParamWidget;
}

class ParamWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit ParamWidget(QWidget *parent = nullptr);
    ~ParamWidget() override;
    std::vector<std::pair<std::string, Vipster::IO::BaseParam>> params;
    void registerParam(const std::string& name,
                       const Vipster::IO::BaseParam& data);
    void clearParams();
    const Vipster::IO::Plugin *curFmt{nullptr};
    Vipster::IO::BaseParam *curParam{nullptr};

private slots:
    void on_paramSel_currentIndexChanged(int index);

    void on_pushButton_clicked();

private:
    Ui::ParamWidget *ui;
    std::map<const Vipster::IO::Plugin*, ParamBase*> formats;
};

#endif // PARAMWIDGET_H
