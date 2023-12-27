#ifndef PSEWIDGET_H
#define PSEWIDGET_H

#include <QWidget>
#include <QListWidgetItem>
#include "../basewidget.h"
#include "vipster/periodictable.h"

namespace Ui {
class PeriodicTableWidget;
}

class PeriodicTableWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PeriodicTableWidget(QWidget *parent = nullptr, bool isGlobal=false);
    ~PeriodicTableWidget();
    void setTable(const Vipster::PeriodicTable* table);

signals:
    void currentElementChanged();

private:
    void setElement(QListWidgetItem* item);

private:
    Ui::PeriodicTableWidget *ui{};
    const bool isGlobal;

    const Vipster::PeriodicTable *table{};
    const std::string* currentName{nullptr};
    const Vipster::Element* currentElement{nullptr};

    template<typename T>
    void registerProperty(QWidget* w, T Vipster::Element::* prop);
};

#endif // PSEWIDGET_H
