#ifndef PYTHONCONSOLE_H
#define PYTHONCONSOLE_H

#include <QPlainTextEdit>
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#pragma pop_macro("slots")

class MainWindow;

class PythonConsole : public QPlainTextEdit
{
    Q_OBJECT

public:
    explicit PythonConsole(QWidget *parent = nullptr);
protected:
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void keyPressEvent(QKeyEvent *e) override;
private:
    std::string getCurCmd();
    size_t curCmd{0};
    std::string tmpCmd{};
    std::vector<std::string> cmdHistory{};
    pybind11::dict locals{};
    int cmdBlock{-1};
};

#endif // PYTHONCONSOLE_H
