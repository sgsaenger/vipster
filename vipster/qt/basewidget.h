#ifndef BASEWIDGET_H
#define BASEWIDGET_H

#include <QWidget>

class MainWindow;

class BaseWidget: public QWidget{
    Q_OBJECT

public:
    BaseWidget(QWidget *parent = nullptr);
    void triggerUpdate(uint8_t change);
    virtual void updateWidget(uint8_t){}
    virtual ~BaseWidget()=default;

protected:
    bool updateTriggered{false};
    MainWindow* master;
};

#endif // BASEWIDGET_H
