#ifndef BASEWIDGET_H
#define BASEWIDGET_H

#include <QWidget>
#include "guiglobals.h"

class MainWindow;

class BaseWidget: public QWidget{
    Q_OBJECT

public:
    BaseWidget(QWidget *parent = nullptr);
    virtual void triggerUpdate(Vipster::GUI::change_t change);
    virtual void updateWidget(Vipster::GUI::change_t){}
    virtual ~BaseWidget()=default;

protected:
    bool updateTriggered{false};
    MainWindow* master;
};

#endif // BASEWIDGET_H
