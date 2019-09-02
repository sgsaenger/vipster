#ifndef PYTHONCONSOLE_H
#define PYTHONCONSOLE_H

#include <QPlainTextEdit>
#pragma push_macro("slots")
#undef slots
#include <pybind11/embed.h>
#pragma pop_macro("slots")
#include "../mainwindow.h"

class PythonConsole : public QPlainTextEdit
{
    Q_OBJECT

public:
    explicit PythonConsole(QWidget *parent = nullptr);
    ~PythonConsole() override;
    void setMaster(MainWindow *master);
protected:
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void keyPressEvent(QKeyEvent *e) override;
private:
    std::string getCurCmd();
    size_t curCmd{0};
    std::string tmpCmd{};
    std::vector<std::string> cmdHistory{};
    QStringList commandbuf{};
    pybind11::scoped_interpreter interp{};
    pybind11::dict locals{};
    int cmdBlock{-1};
    MainWindow* master{nullptr};
};

#endif // PYTHONCONSOLE_H
