#ifndef SCRIPTWIDGET_H
#define SCRIPTWIDGET_H

#include <QWidget>
#include "scripthelp.h"
#include "../mainwindow.h"

namespace Ui {
class ScriptWidget;
}

class ScriptWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit ScriptWidget(QWidget *parent = nullptr);
    ~ScriptWidget() override;
public slots:
    void evalScript();

private slots:
    void on_helpButton_clicked();

private:
    struct OpVec{
        enum class Mode{Direct, Relative, Position, Combination};
        Mode mode{Mode::Direct};
        Vipster::Vec v;
        Vipster::AtomFmt fmt;
        bool m1{false}, m2{false};
        size_t id1, id2;
    };
    friend std::istream& operator>>(std::istream&, std::tuple<OpVec&, bool>);
    struct ScriptOp{
        enum class Mode{None, Rotate, Shift, Mirror, Rename, Select, Define};
        std::string line;
        std::string target;
        Mode mode{Mode::None};
        float f{};
        std::string s1{}, s2{};
        OpVec v1{}, v2{}, v3{};
    };
    std::vector<ScriptOp> parse();
    Vipster::guiChange_t execute(const std::vector<ScriptOp>&,
                                 Vipster::Step&, Vipster::Step::selection&);
    Ui::ScriptWidget *ui;
    ScriptHelp *help;
    static constexpr const char* tooltip=
    "The scripting language is line-based. "
    "The general syntax is:"
    "<p><tt>operator target &lt;arguments&gt;</tt></p>"
    "Only the first three characters of an operator are needed.<br>"
    "The allowed operators are:<ul>"
    "<li><tt><b>sel</b>ect target &lt;selection arguments&gt;</tt>: "
    "change the selection according to arguments (see below)</li>"
    "<li><tt><b>def</b>ine target name &lt;selection arguments&gt;</tt>: "
    "define a named group of atoms according to arguments (see below)</li>"
    "<li><tt><b>shi</b>ft target vec (factor)</tt>: "
    "shift the target by the (optionally scaled) vector</li>"
    "<li><tt><b>rot</b>ate target angle vec (vec)</tt>: "
    "rotate around first vector, optionally shifted by second vector</li>"
    "<li><tt><b>mir</b>ror target vec vec (vec)</tt>: "
    "mirror around a plane given by two vectors, optionally shifted by third vector</li>"
    "<li><tt><b>ren</b>ame target name</tt>: "
    "rename to new atom-type</li>"
    "</ul>"
    "The <tt>target</tt> may be either <tt>all</tt> to match the whole structure, "
    "<tt>sel</tt> to match the current selection at the moment the line is executed, "
    "or an arbitrary name that has been <tt>define</tt>d beforehand.<br><br>"
    "Vectors may be given as:<ul>"
    "<li><tt>(x,y,z)</tt>: Direct vector, evaluated in the <b>set</b> atom format "
    "for the given target</li>"
    "<li><tt>(x,y,z,fmt)</tt>: as above, but explicitely given in a certain format "
    "(one of angstrom, bohr, crystal, alat)</li>"
    "<li><tt>(-)id</tt>: Position vector, the (optionally negated) position "
    "of atom with given id</li>"
    "<li><tt>(-)id1+id2</tt>: Sum vector of two atoms</li>"
    "<li><tt>(-)id1-id2</tt>: Difference vector of two atoms</li>"
    "</ul><b>Attention!</b> Vectors are evaluated for each step, relative formats or "
    "position/difference vectors may not be the same across the trajectory!<br><br>"
    "Selection arguments follow the following BNF-grammar:<p><tt>"
    "Filter ::= (Criterion, {Coupling, Filter}) | (\"(\", Filter, \")\"));<br>"
    "Criterion ::= [\"not \"], (TypeCrit | IdxCrit | PosCrit | CoordCrit);<br>"
    "Coupling ::= [\"!\"], (\"|\" | \"&\" | \"^\");<br>"
    "<br>"
    "TypeCrit ::= \"type \", (Type | TypeList);<br>"
    "TypeList ::= \"[\", Type, {(\" \", Type)}, \"]\";<br>"
    "Type ::= NonWhiteSpace, {NonWhiteSpace};<br>"
    "<br>"
    "IdxCrit ::= \"index \" , (IdxList | IdxRange);<br>"
    "IdxList ::= ( \"[\", IdxRange, {(\" \" IdxRange)}, \"]\");<br>"
    "IdxRange ::= ( Integer, \"-\", Integer) | Integer;<br>"
    "<br>"
    "PosCrit ::= \"pos \", Direction, Format, CompOp, Float;<br>"
    "Direction ::= \"x\" | \"y\" | \"z\";<br>"
    "Format ::= \"a\" | \"b\" | \"c\" | \"d\";<br>"
    "CompOp ::= \">\" | \"<\";<br>"
    "<br>"
    "CoordCrit ::= \"coord \", CompEqOp, Integer;<br>"
    "CompEqOp ::= \"=\" | CompOp;<br>"
    "</tt></p>"
    ;
};

#endif // SCRIPTWIDGET_H
