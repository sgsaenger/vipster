#include "pythonconsole.py.h"
#include <iostream>
#include <QTextBlock>

namespace py = pybind11;
using namespace py::literals;

PythonConsole::PythonConsole(QWidget *parent) :
    QPlainTextEdit(parent)
{
    document()->setDefaultFont(QFont("Courier New", 10));
    setUndoRedoEnabled(false);
    auto ver = py::str(py::module::import("sys").attr("version"));
    textCursor().insertText("Python ");
    textCursor().insertText(std::string{ver}.c_str());
    textCursor().insertText("\nType \"help\" for information about python.\n"
                            "The \"vipster\"-module is preloaded, use \"help"
                            "(vipster)\" for more information.\n>>> ");
}

void PythonConsole::mousePressEvent(QMouseEvent *e)
{

}

void PythonConsole::keyPressEvent(QKeyEvent *e)
{
    switch(e->key()){
    case Qt::Key_Enter:
    case Qt::Key_Return:
        moveCursor(QTextCursor::End);
        if(e->modifiers() & Qt::Modifier::SHIFT){
            textCursor().insertText("\n... ");
        }else{
            textCursor().insertText("\n\n>>> ");
            cmdBlock = document()->lastBlock().blockNumber();
        }
        break;
    case Qt::Key_Up:
    case Qt::Key_Down:
        // scroll through command history or move through multiline command
        // +shift: scroll window?
        break;
    case Qt::Key_Left:
    case Qt::Key_Right:
        // move through current line
        break;
    default:
        moveCursor(QTextCursor::End);
        QPlainTextEdit::keyPressEvent(e);
    }
}


//void PythonConsole::keyReleaseEvent(QKeyEvent *e)
//{
//    QPlainTextEdit::keyReleaseEvent(e);
//}
