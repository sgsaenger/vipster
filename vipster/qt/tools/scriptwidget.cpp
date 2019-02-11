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

std::istream& operator>>(std::istream& is, std::tuple<const Step&, Vec&, bool> dat){
    auto c = static_cast<char>((is >> std::ws).peek());
    const Step& step = std::get<0>(dat);
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
        is >> c >> vec[0] >> c >> vec[1] >> c >> vec[2] >> c;
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
        }
    }else if(c == '-'){
        // negated position vector
        size_t id;
        is >> c >> id;
        vec = -step[id].coord;
    }else{
        // position or difference vector
        size_t id1, id2;
        is >> id1;
        c = static_cast<char>(is.peek());
        if(c == '-'){
            is >> c >> id2;
            vec = step[id1].coord - step[id2].coord;
        }else{
            vec = step[id1].coord;
        }
    }
    return is;
};

void ScriptWidget::evalScript()
{
    guiChange_t change;
    if(ui->trajecCheck->isChecked()){
        change = GuiChange::trajec;
        for(auto& s: master->curMol->getSteps()){
            auto& sel = master->stepdata[&s].sel;
            if(!sel){
                // if step hasn't been loaded before, need to create selection
                sel = std::make_unique<Step::selection>(s.select(SelectionFilter{}));
            }
            change |= evalImpl(s, *sel);
        }
        if(change == GuiChange::trajec){
            change = 0;
        }
    }else{
        change = evalImpl(*master->curStep, *master->curSel);
    }
    triggerUpdate(change);
}

guiChange_t ScriptWidget::evalImpl(Step& step, Step::selection& stepSel)
{
    struct ScriptOp{
        enum class Mode{None, Rotate, Shift, Mirror, Rename, Select, Define};
        std::string target;
        Mode mode{Mode::None};
        float f{};
        std::string s1{}, s2{};
        Vipster::Vec v1{}, v2{}, v3{};
    };
    auto script_str = static_cast<QPlainTextEdit*>(ui->inputEdit)->toPlainText().toStdString();
    auto script = std::stringstream{script_str};
    guiChange_t change{};
    std::map<std::string, Step::selection> definitions{};
    std::vector<ScriptOp> operations{};
    std::string line, op_pre, op(3, ' '), name;
    const bool _false{false}, _true{true};
    while(std::getline(script, line)){
        auto line_stream = std::stringstream{line};
        line_stream >> op_pre;
        std::transform(op_pre.begin(), op_pre.begin()+4, op.begin(), ::tolower);
        // Create op on stack
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
            line_stream >> action.s1;
        }else if(op == "sel"){
            action.mode = ScriptOp::Mode::Select;
            std::getline(line_stream, action.s1);
        }else if(op == "def"){
            action.mode = ScriptOp::Mode::Define;
            line_stream >> action.s1;
            std::getline(line_stream, action.s2);
        }else{
            throw Error("Unknown operator");
        }
    }
    auto execOp = [&definitions, &stepSel](auto& step, const ScriptOp& op){
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
                at.name = op.s1;
            }
            break;
        case ScriptOp::Mode::Select:
            stepSel = step.select(op.s1);
            break;
        case ScriptOp::Mode::Define:
            {
                auto pos = definitions.find(op.s1);
                if(pos == definitions.end()){
                    definitions.insert({op.s1, step.select(op.s2)});
                }else{
                    pos->second = step.select(op.s2);
                }
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
            execOp(stepSel, op);
        }else{
            execOp(definitions.at(op.target), op);
        }
    }
    return change;
}
