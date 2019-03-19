#ifndef PSEWIDGET_H
#define PSEWIDGET_H

#include <QWidget>
#include <QListWidgetItem>
#include "basewidget.h"
#include "pse.h"

namespace Ui {
class PSEWidget;
}

class PSEWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PSEWidget(QWidget *parent = nullptr);
    ~PSEWidget() override;
    void setPSE(Vipster::PseMap* pse);
    void updateWidget(Vipster::guiChange_t) override;

public slots:
    void setEntry(QListWidgetItem* item);

signals:
    void currentEntryChanged();

private:
    Vipster::PseEntry* currentEntry{nullptr};
    Vipster::PseMap* pse{nullptr};
    Ui::PSEWidget *ui;
    bool isGlobal{false};
    template<typename T>
    void registerProperty(QWidget* w, T Vipster::PseEntry::* prop);
};

#endif // PSEWIDGET_H
