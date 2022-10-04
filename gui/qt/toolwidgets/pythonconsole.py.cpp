#include "../mainwindow.h"
#include "pythonconsole.py.h"
#include "glwidget.h"
#include "global.py.h"
#include "vipster/molecule.h"
#include <QTextBlock>
#include <QScrollBar>
#pragma push_macro("slots")
#undef slots
#include <pybind11/embed.h>
#pragma pop_macro("slots")

using namespace Vipster;

PythonConsole::PythonConsole(QWidget *parent) :
    QPlainTextEdit(parent)
{
    // set up console-settings
    document()->setDefaultFont(QFont("Courier New", 10));
    setUndoRedoEnabled(false);
    // print preamble
    auto ver = py::str(py::module::import("sys").attr("version"));
    textCursor().insertText("Python ");
    textCursor().insertText(std::string{ver}.c_str());
    textCursor().insertText("\nType \"help(vipster)\" for more information "
                            "about Vipster-specific functions."
                            "\n>>> ");
    auto vip = py::module::import("vipster");
    // GUI-specific functions
    vip.def("makeScreenshot", [this](const char *fn){master->saveScreenshot(fn);});
    vip.def("viewRotate", [this](float x, float y, float z){
        master->curVP->openGLWidget->rotateViewMat(x, y, z);
    });
    vip.def("viewTranslate", [this](float x, float y, float z){
        master->curVP->openGLWidget->translateViewMat(x, y, z);
    });
    vip.def("viewZoom", [this](float f){
        master->curVP->openGLWidget->zoomViewMat(f);
    });
    // Step-access
    vip.def("curStep", [this](){return master->curStep;}, py::return_value_policy::reference);
    vip.def("getStep", [this](size_t i){
        if(i>=master->curMol->getNstep())
            throw std::range_error("Step-id out of range");
        return master->curMol->getStep(i+1);
    });
    vip.def("setStep", [this](size_t i){
        if(i>=master->curMol->getNstep())
            throw std::range_error("Step-id out of range");
        master->curVP->setStep(i+1);
    });
    vip.def("numStep", [this](){return master->curMol->getNstep();});
    vip.def("curSel", [this](){return master->curSel;}, py::return_value_policy::reference);
    // Molecule-access
    vip.def("curMol", [this](){return master->curMol;}, py::return_value_policy::reference);
    vip.def("getMol", [this](size_t i){
        if(i>=master->molecules.size())
            throw std::range_error("Molecule-id out of range");
        return &*std::next(master->molecules.begin(), i);
    }, py::return_value_policy::reference);
    vip.def("setMol", [this](size_t i){
        if(i>=master->molecules.size())
            throw std::range_error("Molecule-id out of range");
        master->curVP->setMol(i);
    });
    vip.def("numMol", [this](){return master->molecules.size();});
    // Data-access
    vip.def("numData", [this](){return master->data.size();});
    vip.def("getData", [this](size_t i){
        if(i>=master->data.size())
            throw std::range_error("Data-id out of range");
        return std::next(master->data.begin(), i)->get();
    }, py::return_value_policy::reference);
    // import library
    py::exec("import vipster; from vipster import *");
    // overwrite help() to hide interactive mode (breaks console)
    py::exec("def help(*args, **kwds):\n"
             "  if not args and not kwds:\n"
             "    print('Call help(thing) to get information about the python object \"thing\".\\n"
                 "Interactive help is disabled.')\n"
             "  else:\n"
             "    import pydoc\n"
             "    return pydoc.help(*args, **kwds)");
    locals = py::globals();
    cmdBlock = document()->lastBlock().blockNumber();
}

void PythonConsole::setMaster(MainWindow *m)
{
    master = m;
}

PythonConsole::~PythonConsole()
{
}

void PythonConsole::mousePressEvent(QMouseEvent *e)
{
    setReadOnly(true);
    QPlainTextEdit::mousePressEvent(e);
}

void PythonConsole::mouseReleaseEvent(QMouseEvent *e)
{
    QPlainTextEdit::mouseReleaseEvent(e);
    setReadOnly(false);
}

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
    if(e == QKeySequence::Paste){
        // when modifying cmd, make it "active copy"
        curCmd = cmdHistory.size();
        // when we're not in the current cmd block, move there
        if(cursor.blockNumber()<cmdBlock)
            moveCursor(QTextCursor::End);
        int oldCount = blockCount();
        QPlainTextEdit::keyPressEvent(e);
        // expand any newlines that have been inserted
        cursor = textCursor();
        if(blockCount() > oldCount){
            for(int i=0; i<(blockCount()-oldCount); ++i){
                cursor.setPosition(document()->findBlockByNumber(
                                       cursor.blockNumber()-i).position());
                cursor.insertText("... ");
            }
        }
        ensureCursorVisible();
    }else if(e == QKeySequence::Cut){
        if(!cursor.hasSelection()) return;
        if(cursor.selectionStart() <
                document()->findBlockByNumber(cmdBlock).position()){
            if(cursor.selectionEnd() <
                    document()->findBlockByNumber(cmdBlock).position()){
            // if complete selection is outside of current cmd, ignore
                return;
            }else{
            // if selection is partially in current cmd, ignore everything outside
                cursor.setPosition(cursor.selectionEnd());
                cursor.setPosition(document()->findBlockByNumber(cmdBlock).position()+4,
                                   QTextCursor::KeepAnchor);
                setTextCursor(cursor);
            }
        }
        // perform cut
        QPlainTextEdit::keyPressEvent(e);
    }else if(e == QKeySequence::Copy){
        // pass through copy event
        QPlainTextEdit::keyPressEvent(e);
    }else if(e == QKeySequence::DeleteEndOfLine){
        if(cursor.blockNumber()>=cmdBlock){
            if(cursor.selectionStart() >= document()->findBlockByNumber(cmdBlock).position()+4){
                QPlainTextEdit::keyPressEvent(e);
            }
        }
    }else{
        switch(e->key()){
        case Qt::Key_Enter:
        case Qt::Key_Return:
            if(!cursor.atEnd()){
                // insert new line when we're not at end
                cursor.insertText("\n... ");
            }else{
                // otherwise try to run command if complete
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
                        cmdHistory.pop_back();
                        curCmd = cmdHistory.size();
                    }else{
                        // complete command
                        // defer python's stdout
                        auto out = py::module::import("io").attr("StringIO")();
                        auto sys = py::module::import("sys");
                        auto old = sys.attr("stdout");
                        sys.attr("stdout") = out;
                        sys.attr("stderr") = out;
                        // forward compiled code and execute
                        locals["__res__"] = result;
                        py::exec("exec(__res__)", py::globals(), locals);
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
                        // trigger update so GUI knows about changes
                        master->updateWidgets(GUI::stepChanged);
                    }
                } catch (py::error_already_set& e) {
                    cursor.insertText(e.what());
                    cursor.insertText("\n>>> ");
                    cmdBlock = document()->lastBlock().blockNumber();
                }
            }
            break;
        case Qt::Key_Tab:
            // TODO: auto-completion?
            break;
        case Qt::Key_Backspace:
            // when modifying cmd, make it "active copy"
            curCmd = cmdHistory.size();
            if(cursor.blockNumber() >= cmdBlock){
                if(e->modifiers()&Qt::Modifier::CTRL){
                    // do word delete
                    cursor.movePosition(QTextCursor::PreviousWord,
                                        QTextCursor::KeepAnchor);
                    if(cursor.positionInBlock() < 4){
                        // make sure cursor is in valid position after moving a word
                        if(cursor.blockNumber() == cmdBlock){
                            cursor.setPosition(cursor.block().position()+4);
                        }else{
                            cursor.setPosition(cursor.block().position()-1);
                        }
                    }
                    cursor.removeSelectedText();
                }else if(cursor.positionInBlock() > 4){
                    // do regular delete
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
            // when modifying cmd, make it "active copy"
            curCmd = cmdHistory.size();
            if(cursor.blockNumber() >= cmdBlock){
                if(e->modifiers()&Qt::Modifier::CTRL){
                    // do word delete
                    cursor.movePosition(QTextCursor::NextWord,
                                        QTextCursor::KeepAnchor);
                    if(cursor.positionInBlock() < 4){
                        // make sure cursor is in valid position after moving a word
                        cursor.setPosition(cursor.block().position()+4);
                    }
                    cursor.removeSelectedText();
                }else if(cursor.positionInBlock() == cursor.block().length()-1){
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
                    cursor.movePosition(QTextCursor::Up);
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
                    cursor.movePosition(QTextCursor::Down);
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
                auto mode = (e->modifiers()&Qt::Modifier::CTRL) ?
                            QTextCursor::PreviousWord : QTextCursor::Left;
                cursor.movePosition(mode);
                if(cursor.positionInBlock() < 4){
                    if(cursor.blockNumber() == cmdBlock){
                        cursor.setPosition(cursor.block().position()+4);
                    }else{
                        cursor.setPosition(cursor.block().position()-1);
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
                auto mode = (e->modifiers()&Qt::Modifier::CTRL) ?
                            QTextCursor::NextWord : QTextCursor::Right;
                cursor.movePosition(mode);
                if(cursor.positionInBlock() < 4){
                    cursor.setPosition(cursor.block().position()+4);
                }
            }
            break;
        case Qt::Key_Home:
            if(cursor.blockNumber()<cmdBlock){
                cursor.setPosition(document()->findBlockByNumber(cmdBlock).position()+4);
            }else{
                cursor.setPosition(cursor.block().position()+4);
            }
            break;
        case Qt::Key_End:
            if(cursor.blockNumber()<cmdBlock){
                cursor.movePosition(QTextCursor::End);
            }else{
                cursor.movePosition(QTextCursor::EndOfLine);
            }
            break;
        case Qt::Key_PageUp:
            // scrolling without cursor -> early exit
            verticalScrollBar()->triggerAction(QAbstractSlider::SliderPageStepSub);
            return;
        case Qt::Key_PageDown:
            // scrolling without cursor -> early exit
            verticalScrollBar()->triggerAction(QAbstractSlider::SliderPageStepAdd);
            return;
        case Qt::Key_Control:
            // NOP
            return;
        case Qt::Key_Shift:
            // NOP
            return;
        default:
            // when modifying cmd, make it "active copy"
            curCmd = cmdHistory.size();
            // when we're not in the current cmd block, move there
            if(cursor.blockNumber()<cmdBlock)
                moveCursor(QTextCursor::End);
            QPlainTextEdit::keyPressEvent(e);
            ensureCursorVisible();
            return;
        }
        setTextCursor(cursor);
        ensureCursorVisible();
    }
}

