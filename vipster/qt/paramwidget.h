#ifndef PARAMWIDGET_H
#define PARAMWIDGET_H

#include <QWidget>
#include "mainwindow.h"

namespace Ui {
class ParamWidget;
}

class ParamWidget : public QWidget, public Vipster::BaseWidget
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
    Vipster::IO::BaseParam *curParam;
};

#endif // PARAMWIDGET_H
