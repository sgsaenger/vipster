#ifndef LIBVIPSTER_FILTER_H
#define LIBVIPSTER_FILTER_H

#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <set>
#include "atom.h"
#include "bond.h"

/*
 * The Filter mechanism supports multiple criteria
 * as well as combination and nesting
 *
 * Keywords are case insensitive, whereas atom-types are case-sensitive
 *
 * Target grammar (in roughly EBNF):
 *
 * Filter ::= (Criterion, {Coupling, Filter}) | ("(", Filter, ")");
 * Criterion ::= ["not "], (TypeCrit | IdxCrit | PosCrit | CoordCrit);
 * Coupling ::= ["!"], ("|" | "&" | "^");
 *
 * TypeCrit ::= "type ", (Type | TypeList);
 * TypeList ::= "[", Type, {(" ", Type)}, "]";
 * Type ::= NonWhiteSpace, {NonWhiteSpace};
 *
 * IdxCrit ::= "index " , (IdxList | IdxRange);
 * IdxList ::= ( "[", IdxRange, {(" " IdxRange)}, "]");
 * IdxRange ::= ( Integer, "-", Integer) | Integer;
 *
 * PosCrit ::= "pos ", Direction, CompOp, Float;
 * Direction ::= "x" | "y" | "z";
 * CompOp ::= (">" | "<"), ["="];
 *
 * CoordCrit ::= "coord ", CompEqOp, Integer;
 * CompEqOp ::= "=" | CompOp;
 */

namespace Vipster {

// A selected atom is defined via the atom's index and its periodic repetition vector
using SelectionPair = std::pair<size_t, SizeVec>;
using SelectionIndices = std::vector<SelectionPair>;

/* Container for compiled filter, can be converted to/from string
 */
struct SelectionFilter{
    // Filter elements
    enum class Mode:uint8_t{None, Index, Type, Coord, Pos, Group};
    enum Op{NONE=0x0, NOT=0x1,      // first bit negates own op
            PAIR=0x2, NOT_PAIR=0x4, // second bit activates coupling, third bit negates
            AND=0x2, NAND=0x6,
            OR=0xA, NOR=0xE,
            XOR=0x12, XNOR=0x16,
            PAIR_MASK=0x1E,
            UPDATE=0x80};           // marks filter as outdated
    enum Pos{X=0x0, Y=0x1, Z=0x2, DIR_MASK=0x3, // 2  bits for space direction
             P_GT=0x0, P_LT=0x4, P_GEQ=0x8, P_LEQ=0xC,
             P_EQ=0x8, P_CMP_MASK=0xC,  // 2 bit for comp direction
            };
    enum Coord{C_GT=0x0, C_EQ=0x1, C_LT=0x2, C_CMP_MASK=0x3};

    // Interface
    SelectionFilter() = default;
    SelectionFilter(const SelectionFilter& f);
    SelectionFilter(const char *s);
    SelectionFilter(const std::string &s);
    SelectionFilter& operator=(const SelectionFilter& f);
    SelectionFilter& operator=(SelectionFilter&& f);
    operator std::string () const;

    // Data
    Mode mode;
    uint8_t op{Op::UPDATE};
    uint8_t pos;
    uint8_t coord;
    double posVal;
    size_t coordVal;
    SelectionIndices indices;
    std::set<std::string> types;
    std::unique_ptr<SelectionFilter> groupfilter{nullptr};
    std::unique_ptr<SelectionFilter> subfilter{nullptr};
};
std::ostream& operator<<(std::ostream& os, const SelectionFilter& filter);
std::istream& operator>>(std::istream& is, SelectionFilter& filter);

/*
 * Apply filter to a step
 */
template<typename T>
SelectionIndices evalFilter(const T& step, SelectionFilter& filter);

constexpr const char* FilterAbout =
R"--(
<html><head/><body>
<p>A <b><tt>filter</tt></b> is used to pick atoms according to user-defined criteria.</p>
<p>Criteria are as follows:</p>
<ul>
<li><b><tt>type</tt></b>: one or more atom types
<ul>
<li><tt>type C</tt></li>
<li><tt>type [H C N O]</tt></li>
</ul>
</li>
<li><b><tt>index</tt></b>: one or more indices or ranges
<ul>
<li>index 17</li>
<li>index 1-25</li>
<li>index [0 3 5 7-12]</li>
</ul>
</li>
<li><b><tt>pos</tt></b>: relative position (x,y,z to choose axis, &gt;,&gt;=,&lt;,&lt;= to select comparison, followed by target value)
<ul>
<li>pos x&gt;5</li>
<li>pos z &lt; 0.5</li>
</ul>
</li>
<li><b><tt>coord</tt></b>: coordination number
<ul>
<li>coord = 2</li>
<li>coord >0</li>
</ul>
</li>
<li></li>
</ul>
<p>Every criterion can be prefixed with <b><tt>not</tt></b> to invert the selection.</p>
<p>For more complex filters, criteria can be grouped with:</p>
<ul>
<li><b><tt>|</tt></b>: inclusive or, i.e. at least one of</li>
<li><b><tt>!|</tt></b>: inclusive not or, i.e. none of</li>
<li><b><tt>&amp;</tt></b>: and, i.e. both</li>
<li><b><tt>!&amp;</tt></b>: not and, i.e. none or one of, but not both</li>
<li><b><tt>^</tt></b>: exclusive or, i.e. one of, but not both</li>
<li><b><tt>!^</tt></b>: exclusive not or, i.e. either none of or both</li>
</ul>
<p>Furthermore, groupings can be put in parenthesis to simplify logical combinations.</p>
</body></html>)--";
}

#include "filter.tpp"

#endif // LIBVIPSTER_FILTER_H
