#ifndef PARAMWIDGET_H
#define PARAMWIDGET_H

#include <QWidget>
#include "mainwindow.h"

namespace Ui {
class ParamWidget;
}

class ParamBase: public QWidget
{
    Q_OBJECT

public:
    explicit ParamBase(QWidget *parent = nullptr);
    virtual ~ParamBase() = default;
    virtual void setParam(Vipster::BaseParam *p)=0;
};

class ParamWidget : public QWidget, public BaseWidget
{
    Q_OBJECT

public:
    explicit ParamWidget(QWidget *parent = nullptr);
    ~ParamWidget();
    void updateWidget(Vipster::Change change);
    void registerParam(Vipster::IOFmt fmt, const std::string& name);

private slots:
    void on_paramSel_currentIndexChanged(int index);

private:
    Ui::ParamWidget *ui;
    Vipster::BaseParam *curParam;
};

#endif // PARAMWIDGET_H
