#ifndef PSEWIDGET_H
#define PSEWIDGET_H

#include <QWidget>
#include <QListWidgetItem>
#include "mainwindow.h"

namespace Ui {
class PSEWidget;
}

class PSEWidget : public QWidget, public BaseWidget
{
    Q_OBJECT

public:
    explicit PSEWidget(QWidget *parent = nullptr);
    ~PSEWidget();

public slots:
    void setEntry(QListWidgetItem* item);

private:
    Vipster::PseEntry* currentEntry{nullptr};
    Ui::PSEWidget *ui;
    template<typename T>
    void registerProperty(QWidget* w, T Vipster::PseEntry::* prop);
};

#endif // PSEWIDGET_H
