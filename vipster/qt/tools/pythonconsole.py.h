#ifndef PYTHONCONSOLE_H
#define PYTHONCONSOLE_H

#include <QPlainTextEdit>
#pragma push_macro("slots")
#undef slots
#include <pybind11/embed.h>
#pragma pop_macro("slots")

class PythonConsole : public QPlainTextEdit
{
    Q_OBJECT

public:
    explicit PythonConsole(QWidget *parent = nullptr);
    ~PythonConsole() override;
protected:
    void mousePressEvent(QMouseEvent *e) override;
    void keyPressEvent(QKeyEvent *e) override;
private:
    QStringList commandbuf{};
    pybind11::scoped_interpreter *interp;
    pybind11::dict *locals;
    int cmdBlock{-1};
};

#endif // PYTHONCONSOLE_H
