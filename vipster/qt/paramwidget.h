#ifndef PARAMWIDGET_H
#define PARAMWIDGET_H

#include <QWidget>
#include "basewidget.h"
#include "io.h"

namespace Ui {
class ParamWidget;
}

class ParamBase: public QWidget
{
    Q_OBJECT

public:
    explicit ParamBase(QWidget *parent = nullptr);
    virtual ~ParamBase() = default;
    virtual void setParam(Vipster::IO::BaseParam *p)=0;
};

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
};

#endif // PARAMWIDGET_H
