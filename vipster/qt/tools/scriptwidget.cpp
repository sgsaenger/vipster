#include "scriptwidget.h"
#include "ui_scriptwidget.h"
#include "ui_scripthelp.h"
#include "stepsel.h"
#include <QPlainTextEdit>
#include <QMessageBox>

using namespace Vipster;

ScriptWidget::ScriptWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ScriptWidget),
    help(new ScriptHelp)
{
    ui->setupUi(this);
}

ScriptWidget::~ScriptWidget()
{
    delete ui;
    delete help;
}

std::istream& operator>>(std::istream& is, std::tuple<ScriptWidget::OpVec&, bool> dat){
    auto c = static_cast<char>((is >> std::ws).peek());
    auto& vec = std::get<0>(dat);
    bool optional = std::get<1>(dat);
    if(c == '('){
        // explicit vector
        is >> c >> vec.v[0] >> c >> vec.v[1] >> c >> vec.v[2] >> c;
        if(is.fail()){
            throw Error("Explicit vector could not be parsed");
        }
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
            if(fmtmap.find(fmt) == fmtmap.end()){
                throw Error("Unknown vector format: "+fmt);
            }
            vec.mode = ScriptWidget::OpVec::Mode::Relative;
            vec.fmt = fmtmap.at(fmt);
        }else{
            vec.mode = ScriptWidget::OpVec::Mode::Direct;
        }
    }else{
        if(c == '-'){
            vec.m1 = true;
            is >> c;
        }
        is >> vec.id1;
        // if not at end, look ahead
        if(!is.eof()){
            c = static_cast<char>(is.peek());
        }else{
            c = ' ';
        }
        if(c == '-'){
            // difference vector
            vec.m2 = true;
            is >> c >> vec.id2;
            vec.mode = ScriptWidget::OpVec::Mode::Combination;
        }else if(c == '+'){
            // sum vector
            is >> c >> vec.id2;
            vec.mode = ScriptWidget::OpVec::Mode::Combination;
        }else{
            // position vector
            vec.mode = ScriptWidget::OpVec::Mode::Position;
        }
    }
    if(is.fail()){
        if(optional){
            vec.mode = ScriptWidget::OpVec::Mode::Direct;
            vec.v = Vec{};
            return is;
        }
        throw Error("Mandatory vector missing");
    }
    return is;
};

std::istream& operator>>(std::istream& is, std::tuple<float&, bool> dat){
    auto& val = std::get<0>(dat);
    const bool& opt = std::get<1>(dat);
    is >> val;
    if(is.fail() && !opt){
        throw Error("Could not parse mandatory number");
    }
    return is;
}

std::istream& operator>>(std::istream& is, std::tuple<std::string&, bool> dat){
    auto& val = std::get<0>(dat);
    const bool& opt = std::get<1>(dat);
    is >> val;
    if(is.fail() && !opt){
        throw Error("Could not parse mandatory string");
    }
    return is;
}

std::vector<ScriptWidget::ScriptOp> ScriptWidget::parse()
{
    std::vector<ScriptWidget::ScriptOp> operations{};
    auto script_str = static_cast<QPlainTextEdit*>(ui->inputEdit)->toPlainText().toStdString();
    auto script = std::stringstream{script_str};
    std::string line, op_pre, op(3, ' ');
    const bool _false{false}, _true{true};
    try {
        while(std::getline(script, line)){
            auto line_stream = std::stringstream{line};
            // Create op on stack
            operations.push_back({line});
            auto& action = operations.back();
            // parse op
            line_stream >> op_pre;
            std::transform(op_pre.begin(), op_pre.begin()+4, op.begin(), ::tolower);
            // parse arguments
            // TODO: check when and how stream extraction fails
            if(op == "rot"){
                action.mode = ScriptOp::Mode::Rotate;
                line_stream >> action.target;
                line_stream >> std::tie(action.f, _false);
                line_stream >> std::tie(action.v1, _false);
                line_stream >> std::tie(action.v2, _true);
            }else if(op == "shi"){
                action.mode = ScriptOp::Mode::Shift;
                line_stream >> action.target;
                line_stream >> std::tie(action.v1, _false);
                line_stream >> std::tie(action.f, _true);
            }else if(op == "mir"){
                action.mode = ScriptOp::Mode::Mirror;
                line_stream >> action.target;
                line_stream >> std::tie(action.v1, _false);
                line_stream >> std::tie(action.v2, _false);
                line_stream >> std::tie(action.v3, _true);
            }else if(op == "ren"){
                action.mode = ScriptOp::Mode::Rename;
                line_stream >> action.target;
                line_stream >> std::tie(action.s1, _false);
            }else if(op == "sel"){
                action.mode = ScriptOp::Mode::Select;
                std::getline(line_stream, action.s1);
            }else if(op == "def"){
                action.mode = ScriptOp::Mode::Define;
                line_stream >> std::tie(action.s1, _false);
                std::getline(line_stream, action.s2);
            }else{
                throw Error("Unknown operator: "+op);
            }
        }
    } catch (const Error &e) {
        QMessageBox msg{this};
        msg.setText(QString{"Parsing error in script:\n\n"}+line.c_str()+"\n"+e.what());
        msg.exec();
        return {};
    } catch (...) {
        QMessageBox msg{this};
        msg.setText(QString{"Unexpected error when parsing script:\n\n"}+line.c_str());
        msg.exec();
        return {};
    }
    return operations;
}

void ScriptWidget::evalScript()
{
    auto operations = parse();
    guiChange_t change{};
    if(ui->trajecCheck->isChecked()){
        for(auto& s: master->curMol->getSteps()){
            auto& dat = master->stepdata[&s];
            if(!dat.sel){
                // if step hasn't been loaded before, need to create selection
                dat.sel = std::make_unique<Step::selection>(s.select(SelectionFilter{}));
            }
            bool success = execute(operations, s, dat);
            // for current step, save curChange
            if(&s == master->curStep){
                change = curChange;
            }
            // on failure, exit early
            if(!success){
                break;
            }
        }
        if(change) change |= GuiChange::trajec;
    }else{
        execute(operations, *master->curStep, master->stepdata[master->curStep]);
        change = curChange;
    }
    triggerUpdate(change);
}

bool ScriptWidget::execute(const std::vector<ScriptOp>& operations,
                                  Step& step, MainWindow::StepExtras& data)
{
    curChange = guiChange_t{};
    auto mkVec = [&](const OpVec& in)->Vec{
        switch(in.mode){
        case OpVec::Mode::Direct:
            return in.v;
        case OpVec::Mode::Relative:
            return step.formatVec(in.v, in.fmt, step.getFmt());
        case OpVec::Mode::Position:
            if(in.m1){
                return -step.at(in.id1).coord;
            }else{
                return step.at(in.id1).coord;
            }
        case OpVec::Mode::Combination:
        {
            Vec tmp = step.at(in.id1).coord;
            if(in.m1) tmp *= -1;
            if(in.m2){
                tmp -= step.at(in.id2).coord;
            }else{
                tmp += step.at(in.id2).coord;
            }
            return tmp;
        }
        }
    };
    auto execOp = [&](auto& step, const ScriptOp& op){
        switch (op.mode) {
        case ScriptOp::Mode::Rotate:
            step.modRotate(op.f, mkVec(op.v1), mkVec(op.v2));
            curChange |= GuiChange::atoms;
            break;
        case ScriptOp::Mode::Shift:
            step.modShift(mkVec(op.v1), op.f);
            curChange |= GuiChange::atoms;
            break;
        case ScriptOp::Mode::Mirror:
            step.modMirror(mkVec(op.v1), mkVec(op.v2), mkVec(op.v3));
            curChange |= GuiChange::atoms;
            break;
        case ScriptOp::Mode::Rename:
            for(auto& at: step){
                at.name = op.s1;
            }
            curChange |= GuiChange::atoms;
            break;
        case ScriptOp::Mode::Select:
            *data.sel = step.select(op.s1);
            curChange |= GuiChange::selection;
            break;
        case ScriptOp::Mode::Define:
            data.def.insert_or_assign(op.s1, step.select(op.s2));
            curChange |= GuiChange::definitions;
            break;
        default:
            throw Error("Invalid operation");
        }
    };
    for(const auto& op: operations){
        try{
            if(op.target == "all"){
                execOp(step, op);
            }else if(op.target == "sel"){
                execOp(*data.sel, op);
            }else{
                auto def = data.def.find(op.target);
                if(def == data.def.end()){
                    throw Error("Unknown target: "+op.target);
                }
                // make sure that formats match
                def->second.setFmt(step.getFmt());
                execOp(def->second, op);
            }
        } catch (const Error &e) {
            QMessageBox msg{this};
            msg.setText(QString{"Error executing script:\n\n"}+op.line.c_str()+"\n"+e.what());
            msg.exec();
            return false;
        } catch (...) {
            QMessageBox msg{this};
            msg.setText(QString{"Unexpected error when executing line:\n\n"}+op.line.c_str());
            msg.exec();
            return false;
        }
    }
    return true;
}

void ScriptWidget::on_helpButton_clicked()
{
    help->show();
}
