#include "scriptwidget.h"
#include "ui_scriptwidget.h"
#include "stepsel.h"
#include <QPlainTextEdit>

using namespace Vipster;

ScriptWidget::ScriptWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ScriptWidget)
{
    ui->setupUi(this);
}

ScriptWidget::~ScriptWidget()
{
    delete ui;
}

std::istream& operator>>(std::istream& is, std::tuple<Step&, Vec&, bool> dat){
    auto c = static_cast<char>((is >> std::ws).peek());
    Step& step = std::get<0>(dat);
    Vec& vec = std::get<1>(dat);
    bool optional = std::get<2>(dat);
    if(!is.good()){
        if(optional){
            return is;
        }
        throw Error("Mandatory vector missing");
    }
    if(c == '('){
        // explicit vector
        is >> c >> vec[0] >> c >> vec[1] >> c >> vec[2];
        c = static_cast<char>((is >> std::ws).peek());
        if(c != ')'){
            std::string fmt;
            is >> fmt;
            if(fmt.back() == ')'){
                fmt.pop_back();
            }else{
                is >> c;
                if(c != ')'){
                    throw Error("Unterminated vector");
                }
            }
            std::transform(fmt.begin(), fmt.end(), fmt.begin(), ::tolower);
            std::map<std::string, AtomFmt> fmtmap={
                {"angstrom", AtomFmt::Angstrom},
                {"bohr", AtomFmt::Bohr},
                {"crystal", AtomFmt::Crystal},
                {"alat", AtomFmt::Alat}
            };
            vec = step.formatVec(vec, fmtmap.at(fmt), step.getFmt());
        }else{
            is >> c;
        }
    }else if(c == '-'){
        // negated position vector
        size_t id;
        is >> c >> id;
        vec = -step[id].coord;
    }else{
        // position or difference vector
        float id1{}, id2{};
        is >> id1 >> id2;
        if(id2<0){
            vec = step[static_cast<size_t>(id1)].coord - step[static_cast<size_t>(-id2)].coord;
        }else{
            vec = step[static_cast<size_t>(id1)].coord;
        }
    }
    return is;
};

void ScriptWidget::evalScript()
{
    struct ScriptOp{
        enum class Mode{None, Rotate, Shift, Mirror, Rename};
        std::string target;
        Mode mode{Mode::None};
        float f{};
        std::string str{};
        Vipster::Vec v1{}, v2{}, v3{};
    };
    auto script_str = static_cast<QPlainTextEdit*>(ui->inputEdit)->toPlainText().toStdString();
    auto script = std::stringstream{script_str};
    uint8_t change{};
    Step& step = *master->curStep;
    std::map<std::string, Step::selection> definitions{};
    std::vector<ScriptOp> operations{};
    std::string line, op_pre, op(3, ' '), name;
    const bool _false{false}, _true{true};
    while(std::getline(script, line)){
        auto line_stream = std::stringstream{line};
        line_stream >> op_pre;
        std::transform(op_pre.begin(), op_pre.begin()+4, op.begin(), ::tolower);
        if(op == "sel"){
            // Change GUI-selection
            change |= GuiChange::selection;
            std::string sel;
            std::getline(line_stream, sel);
            master->curSel->setFilter(sel);
        }else if(op == "def"){
            // Create group for script-operations
            line_stream >> name;
            std::string sel;
            std::getline(line_stream, sel);
            definitions.insert({name, step.select(sel)});
        }else{
            // Save OPs on stack
            change |= GuiChange::atoms;
            line_stream >> name;
            operations.push_back({name});
            auto& action = operations.back();
            // TODO: check when and how stream extraction fails
            if(op == "rot"){
                action.mode = ScriptOp::Mode::Rotate;
                line_stream >> action.f;
                line_stream >> std::tie(step, action.v1, _false);
                line_stream >> std::tie(step, action.v2, _true);
            }else if(op == "shi"){
                action.mode = ScriptOp::Mode::Shift;
                line_stream >> std::tie(step, action.v1, _false);
                action.f = 1.f;
                line_stream >> action.f;
            }else if(op == "mir"){
                action.mode = ScriptOp::Mode::Mirror;
                line_stream >> std::tie(step, action.v1, _false);
                line_stream >> std::tie(step, action.v2, _false);
                line_stream >> std::tie(step, action.v3, _true);
            }else if(op == "ren"){
                action.mode = ScriptOp::Mode::Rename;
                line_stream >> action.str;
            }else{
                throw Error("Unknown operator");
            }
        }
    }
    auto execOp = [](auto& step, const ScriptOp& op){
        switch (op.mode) {
        case ScriptOp::Mode::Rotate:
            step.modRotate(op.f, op.v1, op.v2);
            break;
        case ScriptOp::Mode::Shift:
            step.modShift(op.v1, op.f);
            break;
        case ScriptOp::Mode::Mirror:
            step.modMirror(op.v1, op.v2, op.v3);
            break;
        case ScriptOp::Mode::Rename:
            for(auto& at: step){
                at.name = op.str;
            }
            break;
        default:
            throw Error("Invalid operation");
        }
    };
    for(const auto& op: operations){
        if(op.target == "all"){
            execOp(step, op);
        }else if(op.target == "sel"){
            execOp(*master->curSel, op);
        }else{
            execOp(definitions.at(op.target), op);
        }
    }
    triggerUpdate(change);
}
