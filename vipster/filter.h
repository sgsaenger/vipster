#ifndef LIBVIPSTER_FILTER_H
#define LIBVIPSTER_FILTER_H

#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <set>
#include "atom.h"
#include "bond.h"

namespace Vipster {

using SelectionPair = std::pair<size_t, SizeVec>;
using SelectionIndices = std::vector<SelectionPair>;

/*
 * Wrap multiple filter criteria without polymorphism
 * also represents hierarchically coupled filter chains
 * Target grammar:
 * Keywords are case insensitive, whereas atom-types are case-sensitive
 *
 * Filter ::= (Criterion, {Coupling, Filter}) | ("(", Filter, ")"));
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
 * PosCrit ::= "pos ", Direction, Format, CompOp, Float;
 * Direction ::= "x" | "y" | "z";
 * CompOp ::= (">" | "<"), ["="];
 *
 * CoordCrit ::= "coord ", CompEqOp, Integer;
 * CompEqOp ::= "=" | CompOp;
 */
struct SelectionFilter;
std::ostream& operator<<(std::ostream& os, const SelectionFilter& filter);
std::istream& operator>>(std::istream& is, SelectionFilter& filter);

struct SelectionFilter{
    SelectionFilter() = default;
    SelectionFilter(const SelectionFilter& f) {
        *this = f;
    }
    SelectionFilter(const char *s){
        std::stringstream ss{s};
        ss >> *this;
    }
    SelectionFilter(const std::string &s){
        std::stringstream ss{s};
        ss >> *this;
    }
    SelectionFilter& operator=(const SelectionFilter& f){
        mode = f.mode; op = f.op | Op::UPDATE;
        pos = f.pos; posVal = f.posVal;
        coord = f.coord; coordVal = f.coordVal;
        indices = f.indices; types = f.types;
        if(f.subfilter){
            subfilter = std::make_unique<SelectionFilter>(*f.subfilter);
        }
        if(f.groupfilter){
            groupfilter = std::make_unique<SelectionFilter>(*f.groupfilter);
        }
        return *this;
    }
    SelectionFilter& operator=(SelectionFilter&& f){
        mode = f.mode;
        op = f.op | Op::UPDATE;
        pos = f.pos;
        posVal = f.posVal;
        coord = f.coord;
        coordVal = f.coordVal;
        indices = std::move(f.indices);
        types = std::move(f.types);
        subfilter = std::move(f.subfilter);
        groupfilter = std::move(f.groupfilter);
        return *this;
    }
    operator std::string () const
    {
        std::stringstream ss{};
        ss << *this;
        return ss.str();
    }
    enum class Mode:uint8_t{None, Index, Type, Coord, Pos, Group};
    enum Op{NONE=0x0, NOT=0x1, // first bit negates own op
            PAIR=0x2, NOT_PAIR=0x4, // second bit activates coupling, third bit negates
            AND=0x2, NAND=0x6,
            OR=0xA, NOR=0xE,
            XOR=0x12, XNOR=0x16,
            PAIR_MASK=0x1E,
            UPDATE=0x80};
    enum Pos{X=0x0, Y=0x1, Z=0x2, DIR_MASK=0x3, // 2  bits for space direction
             P_GT=0x0, P_LT=0x4, P_GEQ=0x8, P_LEQ=0xC,
             P_EQ=0x8, P_CMP_MASK=0xC,  // 2 bit for comp direction
            };
    enum Coord{C_GT=0x0, C_EQ=0x1, C_LT=0x2, C_CMP_MASK=0x3};
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

/*
 * recursively evaluate the filters
 */
template<typename T>
static SelectionIndices evalType(const T& step, const SelectionFilter& filter){
    SelectionIndices tmp;
    size_t idx{0};
    for(const auto& at: step){
        for(const auto& type: filter.types){
            if(at.name == type){
                tmp.emplace_back(idx, SizeVec{});
                break;
            }
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static SelectionIndices evalIdx(const T& step, const SelectionFilter& filter){
    SelectionIndices tmp;
    auto nat = step.getNat();
    for(const auto& p: filter.indices){
        if(p.first<nat){
            tmp.push_back(p);
        }
    }
    return tmp;
}

template<typename T>
static SelectionIndices evalPos(const T& step, const SelectionFilter& filter){
    SelectionIndices tmp;
    std::size_t idx{0};
    auto cmp = [&filter](const Vec& at){
        size_t dir = (filter.pos & filter.DIR_MASK);
        switch(filter.pos & filter.P_CMP_MASK){
        case SelectionFilter::P_GT:
            return at[dir] > filter.posVal;
        case SelectionFilter::P_LT:
            return at[dir] < filter.posVal;
        case SelectionFilter::P_GEQ:
            return at[dir] >= filter.posVal;
        case SelectionFilter::P_LEQ:
            return at[dir] <= filter.posVal;
        }
    };
    for(const auto& at: step){
        if(cmp(at.coord)){
            tmp.emplace_back(idx, SizeVec{});
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static SelectionIndices evalCoord(const T& step, const SelectionFilter& filter){
    SelectionIndices tmp;
    // TODO: move functionality to step?
    std::vector<size_t> coord_numbers(step.getNat());
    for(const Bond& b: step.getBonds()){
        coord_numbers[b.at1] += 1;
        coord_numbers[b.at2] += 1;
    }
    auto cmp = [&filter](const size_t c){
        auto cmp_op = filter.coord & filter.C_CMP_MASK;
        if(cmp_op == filter.C_GT){
            return c > filter.coordVal;
        }else if(cmp_op == filter.C_EQ){
            return c == filter.coordVal;
        }else{
            return c < filter.coordVal;
        }
    };
    for(size_t i=0; i<step.getNat(); ++i){
        if(cmp(coord_numbers[i])){
            tmp.emplace_back(i, SizeVec{});
        }
    }
    return tmp;
}

template<typename T>
static SelectionIndices invertSel(const T& step, const SelectionIndices& in){
    SelectionIndices out;
    const auto nat = step.getNat();
    out.reserve(nat);
    for(size_t i=0; i<nat; ++i){
        if(std::find_if(in.begin(), in.end(),
                        [&i](const auto& pair){
                            return pair.first == i;
                        }) == in.end()){
            out.emplace_back(i, SizeVec{});
        }
    }
    return out;
}

template<typename T>
static SelectionIndices evalSubFilter(const T& step,
                                  const SelectionFilter& filter,
                                  SelectionFilter& subfilter,
                                  const SelectionIndices& parent){
    using Op = SelectionFilter::Op;
    SelectionIndices child = evalFilter(step, subfilter);
    SelectionIndices tmp(parent.size()+child.size());
    SelectionIndices::iterator it;
    if((filter.op & Op::XOR) == Op::XOR){
        it = std::set_symmetric_difference(parent.begin(), parent.end(),
                                           child.begin(), child.end(),
                                           tmp.begin());
    }else if((filter.op & Op::OR) == Op::OR){
        it = std::set_union(parent.begin(), parent.end(),
                            child.begin(), child.end(),
                            tmp.begin());
    }else if((filter.op & Op::AND) == Op::AND){
        it = std::set_intersection(parent.begin(), parent.end(),
                                   child.begin(), child.end(),
                                   tmp.begin());
    }else{
        throw Error("Unknown coupling operator "+std::to_string(filter.op & Op::PAIR_MASK));
    }
    tmp.resize(static_cast<size_t>(it - tmp.begin()));
    if(filter.op & Op::NOT_PAIR){
        return invertSel(step, tmp);
    }
    return tmp;
}

template<typename T>
SelectionIndices evalFilter(const T& step, SelectionFilter& filter)
{
    SelectionIndices tmp;
    switch(filter.mode){
    case SelectionFilter::Mode::Group:
        tmp = evalFilter(step, *filter.groupfilter);
        break;
    case SelectionFilter::Mode::Type:
        tmp = evalType(step, filter);
        break;
    case SelectionFilter::Mode::Index:
        tmp = evalIdx(step, filter);
        break;
    case SelectionFilter::Mode::Pos:
        tmp = evalPos(step, filter);
        break;
    case SelectionFilter::Mode::Coord:
        tmp = evalCoord(step, filter);
        break;
    default:
        return tmp;
    }
    if(filter.op & filter.NOT){
        tmp = invertSel(step, tmp);
    }
    if(filter.op & filter.PAIR){
        tmp = evalSubFilter(step, filter, *filter.subfilter, tmp);
    }
    filter.op &= ~filter.UPDATE;
    return tmp;
}

constexpr const char* FilterAbout =
        "<html><head/><body>"
        "<p>A <b><tt>filter</tt></b> is used to pick atoms according to user-defined criteria.</p>"
        "<p>Criteria are as follows:</p>"
        "<ul>"
        "<li><b><tt>type</tt></b>: one or more atom types"
        "<ul>"
        "<li><tt>type C</tt></li>"
        "<li><tt>type [H C N O]</tt></li>"
        "</ul>"
        "</li>"
        "<li><b><tt>index</tt></b>: one or more indices or ranges"
        "<ul>"
        "<li>index 17</li>"
        "<li>index 1-25</li>"
        "<li>index [0 3 5 7-12]</li>"
        "</ul>"
        "</li>"
        "<li><b><tt>pos</tt></b>: relative position (x,y,z to choose axis, &gt;,&gt;=,&lt;,&lt;= to select comparison, followed by target value)"
        "<ul>"
        "<li>pos xa&gt;5</li>"
        "<li>pos z c &lt; 0.5</li>"
        "</ul>"
        "</li>"
        "<li><b><tt>coord</tt></b>: coordination number"
        "<ul>"
        "<li>coord = 2</li>"
        "<li>coord >0</li>"
        "</ul>"
        "</li>"
        "<li></li>"
        "</ul>"
        "<p>Every criterion can be prefixed with <b><tt>not</tt></b> to invert the selection.</p>"
        "<p>For more complex filters, criteria can be grouped with:</p>"
        "<ul>"
        "<li><b><tt>|</tt></b>: inclusive or, i.e. at least one of</li>"
        "<li><b><tt>!|</tt></b>: inclusive not or, i.e. none of</li>"
        "<li><b><tt>&amp;</tt></b>: and, i.e. both</li>"
        "<li><b><tt>!&amp;</tt></b>: not and, i.e. none or one of, but not both</li>"
        "<li><b><tt>^</tt></b>: exclusive or, i.e. one of, but not both</li>"
        "<li><b><tt>!^</tt></b>: exclusive not or, i.e. either none of or both</li>"
        "</ul>"
        "<p>Furthermore, groupings can be put in parenthesis to simplify logical combinations.</p>"
        "</body></html>";
}

#endif // LIBVIPSTER_FILTER_H
