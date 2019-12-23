#ifndef PARAMWIDGET_H
#define PARAMWIDGET_H

#include "../basewidget.h"
#include "fileio.h"

#include <QTreeWidget>

namespace Ui {
class ParamWidget;
}

class ParamWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit ParamWidget(QWidget *parent = nullptr);
    ~ParamWidget() override;
    std::vector<std::pair<std::string, Vipster::IO::Parameter>> params;
    void registerParam(const std::string& name,
                       const Vipster::IO::Parameter& data);
    void clearParams();
    Vipster::IO::Parameter *curParam{nullptr};

private slots:
    void on_paramSel_currentIndexChanged(int index);

    void on_pushButton_clicked();

private:
    Ui::ParamWidget *ui;
    // save state for text edit stack
    std::string curVec{};
    // save state for treeview
    std::string curKey{};
    QTreeWidgetItem *curItem{nullptr};
    std::map<std::string, std::string>* curNL{nullptr};
    void setupText(const QVector<QString> &vectors);
    QTreeWidget* setupTree();
};

#endif // PARAMWIDGET_H
