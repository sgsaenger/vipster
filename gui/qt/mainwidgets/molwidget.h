#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QScrollArea>

namespace Ui {
class MolWidget;
}

class MolWidget : public QScrollArea
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = nullptr);
    ~MolWidget() override;

private:
    Ui::MolWidget *ui;
};

#endif // MOLWIDGET_H
