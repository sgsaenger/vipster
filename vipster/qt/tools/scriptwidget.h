#ifndef SCRIPTWIDGET_H
#define SCRIPTWIDGET_H

#include <QWidget>
#include "../mainwindow.h"

namespace Ui {
class ScriptWidget;
}

class ScriptWidget : public QWidget, public BaseWidget
{
    Q_OBJECT

public:
    explicit ScriptWidget(QWidget *parent = nullptr);
    ~ScriptWidget() override;
public slots:
    void evalScript();

private:
    Ui::ScriptWidget *ui;
};

#endif // SCRIPTWIDGET_H
