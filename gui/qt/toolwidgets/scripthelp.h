#ifndef SCRIPTHELP_H
#define SCRIPTHELP_H

#include <QDialog>

namespace Ui {
class ScriptHelp;
}

class ScriptHelp : public QDialog
{
    Q_OBJECT

public:
    explicit ScriptHelp(QWidget *parent = nullptr);
    ~ScriptHelp();

private:
    Ui::ScriptHelp *ui;
};

#endif // SCRIPTHELP_H
