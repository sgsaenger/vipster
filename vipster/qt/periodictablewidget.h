#ifndef PSEWIDGET_H
#define PSEWIDGET_H

#include <QWidget>
#include <QListWidgetItem>
#include "basewidget.h"
#include "periodictable.h"

namespace Ui {
class PeriodicTableWidget;
}

class PeriodicTableWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PeriodicTableWidget(QWidget *parent = nullptr);
    ~PeriodicTableWidget() override;
    void setTable(Vipster::PeriodicTable* table);
    void updateWidget(Vipster::GUI::change_t) override;

public slots:
    void setEntry(QListWidgetItem* item);
    void changeEntry();

signals:
    void currentEntryChanged();

private slots:
    void on_helpButton_clicked();

private:
    const std::string* currentName{nullptr};
    Vipster::Element* currentElement{nullptr};
    Vipster::PeriodicTable* table{nullptr};
    Ui::PeriodicTableWidget *ui;
    bool isGlobal{false};
    template<typename T>
    void registerProperty(QWidget* w, T Vipster::Element::* prop);
};

#endif // PSEWIDGET_H
