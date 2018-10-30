#ifndef BASEWIDGET_H
#define BASEWIDGET_H

#include <QWidget>
#include "../common/guiglobals.h"

class MainWindow;

class BaseWidget: public QWidget{
    Q_OBJECT

public:
    BaseWidget(QWidget *parent = nullptr);
    void triggerUpdate(Vipster::guiChange_t change);
    virtual void updateWidget(Vipster::guiChange_t){}
    virtual ~BaseWidget()=default;

protected:
    bool updateTriggered{false};
    MainWindow* master;
};

#endif // BASEWIDGET_H
