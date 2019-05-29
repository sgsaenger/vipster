#ifndef PYTHONCONSOLE_H
#define PYTHONCONSOLE_H

#include <QPlainTextEdit>
#pragma push_macro("slots")
#undef slots
#include "pybind11/embed.h"
#pragma pop_macro("slots")

class PythonConsole : public QPlainTextEdit
{
    Q_OBJECT

public:
    explicit PythonConsole(QWidget *parent = nullptr);
protected:
    void mousePressEvent(QMouseEvent *e) override;
    void keyPressEvent(QKeyEvent *e) override;
//    void keyReleaseEvent(QKeyEvent *e) override;
private:
    QStringList commandbuf{};
    int cmdBlock{-1};
    pybind11::scoped_interpreter interp{};
};

#endif // PYTHONCONSOLE_H
