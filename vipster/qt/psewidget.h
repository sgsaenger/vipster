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
    void setPSE(Vipster::PseMap* pse);
    void updateWidget(uint8_t) override;

public slots:
    void setEntry(QListWidgetItem* item);

signals:
    void currentEntryChanged();

private:
    Vipster::PseEntry* currentEntry{nullptr};
    Vipster::PseMap* pse{nullptr};
    Ui::PSEWidget *ui;
    template<typename T>
    void registerProperty(QWidget* w, T Vipster::PseEntry::* prop);
};

#endif // PSEWIDGET_H
