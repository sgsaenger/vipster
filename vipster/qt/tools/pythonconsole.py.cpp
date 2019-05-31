#include "pythonconsole.py.h"
#include <iostream>
#include <pybind11/iostream.h>
#include <QTextBlock>

namespace py = pybind11;
using namespace py::literals;

PythonConsole::PythonConsole(QWidget *parent) :
    QPlainTextEdit(parent)
{
    // set up console-settings
    document()->setDefaultFont(QFont("Courier New", 10));
    setUndoRedoEnabled(false);
    // create static python objects
    interp = new py::scoped_interpreter{};
    locals = new py::dict{};
    // print preamble
    auto ver = py::str(py::module::import("sys").attr("version"));
    textCursor().insertText("Python ");
    textCursor().insertText(std::string{ver}.c_str());
    textCursor().insertText("\nType \"help\" for information about python."
                            "\nType \"help(vipster)\" for more information "
                            "about Vipster-specific functions."
                            "\n>>> ");
    cmdBlock = document()->lastBlock().blockNumber();
}

PythonConsole::~PythonConsole()
{
    delete locals;
    delete interp;
}

void PythonConsole::mousePressEvent(QMouseEvent *e)
{

}

void PythonConsole::keyPressEvent(QKeyEvent *e)
{
    auto cursor = textCursor();
    switch(e->key()){
    case Qt::Key_Enter:
    case Qt::Key_Return:
        moveCursor(QTextCursor::End);
        {
            std::string cmd =
                document()->findBlockByNumber(cmdBlock).
                    text().toStdString().c_str()+4;
            for(int i=cmdBlock+1; i<document()->blockCount(); ++i){
                const auto& block = document()->findBlockByNumber(i);
                // copy block without first 4 chars
                cmd += '\n';
                cmd += block.text().toStdString().c_str()+4;
            }
            cursor.insertText("\n");
            try {
                auto code = py::module::import("code");
                auto result = code.attr("compile_command")(cmd);
                if(result.is_none()){
                    // incomplete command
                    cursor.insertText("... ");
                }else{
                    // complete command
                    auto out = py::module::import("io").attr("StringIO")();
                    auto sys = py::module::import("sys");
                    auto old = sys.attr("stdout");
                    sys.attr("stdout") = out;
                    (*locals)["__res__"] = result;
                    py::exec("exec(__res__)", py::globals(), *locals);
                    auto outval = out.attr("getvalue")();
                    if(py::len(outval)){
                        std::string output = py::str(outval);
                        cursor.insertText(output.c_str());
                    }
                    cursor.insertText("\n>>> ");
                    cmdBlock = document()->lastBlock().blockNumber();
                    sys.attr("stdout") = old;
                }
            } catch (py::error_already_set& e) {
                cursor.insertText(e.what());
                cursor.insertText("\n>>> ");
                cmdBlock = document()->lastBlock().blockNumber();
            }
            ensureCursorVisible();
        }
        break;
    case Qt::Key_Tab:
        // TODO: auto-completion?
        break;
    case Qt::Key_Backspace:
        // TODO: disallow deleting selection
        if(cursor.block().blockNumber() >= cmdBlock){
            if(cursor.positionInBlock() > 4){
                // in cmd, do regular delete
                QPlainTextEdit::keyPressEvent(e);
            }else if(cursor.block().blockNumber() > cmdBlock){
                // in continuation line, remove
                cursor.deletePreviousChar();
                cursor.deletePreviousChar();
                cursor.deletePreviousChar();
                cursor.deletePreviousChar();
                cursor.deletePreviousChar();
            }
            // else we are in first cmd line, do nothing
        }
        break;
    case Qt::Key_Delete:
        // TODO
        // delete only until last >>>
        std::cout << cursor.position() << std::endl;
        std::cout << cursor.positionInBlock() << std::endl;
        std::cout << cursor.block().blockNumber() << std::endl;
        std::cout << cmdBlock << std::endl;
        break;
    case Qt::Key_Up:
    case Qt::Key_Down:
        // TODO
        // scroll through command history or move through multiline command
        // +shift: scroll window?
        break;
    case Qt::Key_Left:
    case Qt::Key_Right:
        // TODO
        // move through current line
        break;
    case Qt::Key_Home:
        // TODO
        cursor.setPosition(cursor.block().position()+4);
        break;
    case Qt::Key_End:
        // TODO
        break;
    case Qt::Key_Insert:
        // TODO
        break;
    default:
        moveCursor(QTextCursor::End);
        QPlainTextEdit::keyPressEvent(e);
    }
    // TODO: limit keycombos like ctrl+a?
}

