#include "pythonconsole.py.h"
#include <iostream>
#include <pybind11/iostream.h>
#include <QTextBlock>

namespace py = pybind11;
using namespace py::literals;

PythonConsole::PythonConsole(QWidget *parent) :
    QPlainTextEdit(parent),
    interp(new py::scoped_interpreter{}),
    locals(new py::dict{})
{
    // set up console-settings
    document()->setDefaultFont(QFont("Courier New", 10));
    setUndoRedoEnabled(false);
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

//void PythonConsole::mousePressEvent(QMouseEvent *e)
//{

//}

std::string PythonConsole::getCurCmd()
{
    if(curCmd >= cmdHistory.size()){
        std::string cmd =
            document()->findBlockByNumber(cmdBlock).
                text().toStdString().c_str()+4;
        for(int i=cmdBlock+1; i<blockCount(); ++i){
            const auto& block = document()->findBlockByNumber(i);
            // copy block without first 4 chars
            cmd += '\n';
            cmd += block.text().toStdString().c_str()+4;
        }
        return cmd;
    }else{
        return cmdHistory[curCmd];
    }
}

void PythonConsole::keyPressEvent(QKeyEvent *e)
{
    auto cursor = textCursor();
    auto move = e->modifiers()&Qt::Modifier::SHIFT ?
                QTextCursor::KeepAnchor : QTextCursor::MoveAnchor;
    switch(e->key()){
    case Qt::Key_Enter:
    case Qt::Key_Return:
        moveCursor(QTextCursor::End);
        cmdHistory.push_back(getCurCmd());
        curCmd = cmdHistory.size();
        tmpCmd.clear();
        cursor.insertText("\n");
        try {
            auto code = py::module::import("code");
            auto result = code.attr("compile_command")(cmdHistory.back());
            if(result.is_none()){
                // incomplete command
                cursor.insertText("... ");
            }else{
                // complete command
                // defer python's stdout
                auto out = py::module::import("io").attr("StringIO")();
                auto sys = py::module::import("sys");
                auto old = sys.attr("stdout");
                sys.attr("stdout") = out;
                // forward compiled code and execute
                (*locals)["__res__"] = result;
                py::exec("exec(__res__)", py::globals(), *locals);
                // handle output
                auto outval = out.attr("getvalue")();
                if(py::len(outval)){
                    std::string output = py::str(outval);
                    cursor.insertText(output.c_str());
                }
                cursor.insertText("\n>>> ");
                // reset state for new execution
                cmdBlock = document()->lastBlock().blockNumber();
                sys.attr("stdout") = old;
            }
        } catch (py::error_already_set& e) {
            cursor.insertText(e.what());
            cursor.insertText("\n>>> ");
            cmdBlock = document()->lastBlock().blockNumber();
        }
        ensureCursorVisible();
        break;
    case Qt::Key_Tab:
        // TODO: auto-completion?
        break;
    case Qt::Key_Backspace:
        if(cursor.blockNumber() >= cmdBlock){
            if(cursor.positionInBlock() > 4){
                // in cmd, do regular delete
                cursor.deletePreviousChar();
            }else if(cursor.block().blockNumber() > cmdBlock){
                // in continuation line, remove until end of last line
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
        if(cursor.blockNumber() >= cmdBlock){
            if(cursor.positionInBlock() == cursor.block().length()-1){
                // if continuation line exists, delete prefix
                if(cursor.blockNumber() < blockCount()){
                    cursor.deleteChar();
                    cursor.deleteChar();
                    cursor.deleteChar();
                    cursor.deleteChar();
                    cursor.deleteChar();
                }
            }else{
                // regular delete
                cursor.deleteChar();
            }
        }
        break;
    case Qt::Key_Up:
        if(cursor.blockNumber() < cmdBlock){
            // move to current cmd-block
            cursor.movePosition(QTextCursor::End);
        }else{
            if(cursor.blockNumber() > cmdBlock){
                // move inside of block
                cursor.movePosition(QTextCursor::Up, move);
            }else{
                // fetch previous from history
                if(curCmd > 0){
                    if(curCmd >= cmdHistory.size()){
                        // save current cmd if edited
                        tmpCmd = getCurCmd();
                        curCmd = cmdHistory.size()-1;
                    }else{
                        curCmd--;
                    }
                    cursor.setPosition(document()->findBlockByNumber(cmdBlock).position()+4);
                    cursor.movePosition(QTextCursor::End, QTextCursor::KeepAnchor);
                    cursor.removeSelectedText();
                    QString tmp{};
                    for(const auto& c: cmdHistory[curCmd]){
                        tmp += c;
                        if(c == '\n'){
                            tmp += "... ";
                        }
                    }
                    cursor.insertText(tmp);
                }
            }
        }
        break;
    case Qt::Key_Down:
        if(cursor.blockNumber() < cmdBlock){
            // move to current cmd-block
            cursor.movePosition(QTextCursor::End);
        }else{
            if(cursor.blockNumber() < blockCount()-1){
                // move inside of block
                cursor.movePosition(QTextCursor::Down, move);
            }else{
                // fetch later from history
                if(curCmd < cmdHistory.size()){
                    curCmd++;
                    cursor.setPosition(document()->findBlockByNumber(cmdBlock).position()+4);
                    cursor.movePosition(QTextCursor::End, QTextCursor::KeepAnchor);
                    cursor.removeSelectedText();
                    QString tmp{};
                    const std::string &source = curCmd < cmdHistory.size() ?
                                cmdHistory[curCmd] : tmpCmd;
                    for(const auto& c: source){
                        tmp += c;
                        if(c == '\n'){
                            tmp += "... ";
                        }
                    }
                    cursor.insertText(tmp);
                }
            }
        }
        break;
    case Qt::Key_Left:
        if(cursor.blockNumber() < cmdBlock){
            // move to current cmd-block
            cursor.movePosition(QTextCursor::End);
        }else{
            // move through current cmd-block
            auto mode = e->modifiers()&Qt::Modifier::CTRL ?
                        QTextCursor::PreviousWord : QTextCursor::Left;
            cursor.movePosition(mode, move);
            if(cursor.positionInBlock() < 4){
                if(cursor.blockNumber() == cmdBlock){
                    cursor.setPosition(cursor.block().position()+4, move);
                }else{
                    cursor.setPosition(cursor.block().position()-1, move);
                }
            }
        }
        break;
    case Qt::Key_Right:
        if(cursor.blockNumber() < cmdBlock){
            // move to current cmd-block
            cursor.setPosition(document()->findBlockByNumber(cmdBlock).position()+4);
        }else{
            // move through current cmd-block
            auto mode = e->modifiers()&Qt::Modifier::CTRL ?
                        QTextCursor::NextWord : QTextCursor::Right;
            cursor.movePosition(mode, move);
            if(cursor.positionInBlock() < 4){
                cursor.setPosition(cursor.block().position()+4, move);
            }
        }
        break;
    case Qt::Key_Home:
        if(cursor.blockNumber()<cmdBlock){
            cursor.setPosition(document()->findBlockByNumber(cmdBlock).position()+4);
        }else{
            cursor.setPosition(cursor.block().position()+4, move);
        }
        break;
    case Qt::Key_End:
        if(cursor.blockNumber()<cmdBlock){
            cursor.movePosition(QTextCursor::End);
        }else{
            cursor.movePosition(QTextCursor::EndOfLine, move);
        }
        break;
    case Qt::Key_PageUp:
    case Qt::Key_PageDown:
        //TODO
        break;
    default:
        // when modifying cmd, make it "active copy"
        std::cout << curCmd << std::endl;
        curCmd = cmdHistory.size();
        std::cout << curCmd << '\n' << std::endl;
        // when we're not in the current cmd block, move there
        if(cursor.blockNumber()<cmdBlock)
            moveCursor(QTextCursor::End);
        QPlainTextEdit::keyPressEvent(e);
        return;
    }
    setTextCursor(cursor);
    // TODO: limit keycombos like ctrl+a?
    // TODO: modifiers for cursor movement
    // TODO: handle selection somehow
}

