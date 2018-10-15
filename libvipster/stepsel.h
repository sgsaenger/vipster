#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <set>
#include "atom.h"
#include "bond.h"
#include "cell.h"

namespace Vipster {
// Forward declarations
template<typename T>
class StepConst;
class Step;

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
 * Format ::= "a" | "b" | "c" | "d";
 * CompOp ::= ">" | "<";
 *
 * CoordCrit ::= "coord ", CompEqOp, Integer;
 * CompEqOp ::= "=" | CompOp;
 */
struct SelectionFilter;
std::ostream& operator<<(std::ostream& os, const SelectionFilter& filter);
std::istream& operator>>(std::istream& is, SelectionFilter& filter);
//TODO: needs to be templated, maybe member-function. ATM copy-creates Step EVERY time
//std::vector<size_t> evalFilter(const Step& step, SelectionFilter& filter);

struct SelectionFilter{
    SelectionFilter() = default;
    SelectionFilter(const SelectionFilter& f) {
        *this = f;
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
        mode = f.mode; op = f.op | Op::UPDATE;
        pos = f.pos; posVal = f.posVal;
        coord = f.coord; coordVal = f.coordVal;
        indices = f.indices; types = f.types;
        subfilter = std::move(f.subfilter);
        groupfilter = std::move(f.groupfilter);
        return *this;
    }
    std::string toStr() const
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
    enum Pos{BOHR=0x0, ANG=0x1, CRYS=0x2, CDM=0x3, FMT_MASK=0x3,// 2 bits for format
             P_GT=0x0, P_LT=0x4, P_CMP_MASK=0x4,  // 1 bit for comp direction
             X=0x0, Y=0x8, Z=0x10, DIR_MASK=0x18, // 2  bits for space direction
            };
    enum Coord{C_GT=0x0, C_EQ=0x1, C_LT=0x2, C_CMP_MASK=0x3};
    Mode mode;
    uint8_t op{Op::UPDATE};
    uint8_t pos;
    uint8_t coord;
    float posVal;
    size_t coordVal;
    std::set<size_t> indices;
    std::set<std::string> types;
    std::unique_ptr<SelectionFilter> groupfilter{nullptr};
    std::unique_ptr<SelectionFilter> subfilter{nullptr};
};

/*
 * recursively evaluate the filters
 */
template<typename T>
static std::vector<size_t> evalType(const T& step, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    size_t idx{0};
    for(const auto& at: step){
        for(const auto& type: filter.types){
            if(at.name == type){
                tmp.push_back(idx);
                break;
            }
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static std::vector<size_t> evalIdx(const T& step, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    auto nat = step.getNat();
    for(const auto& i:filter.indices){
        if(i<nat){
            tmp.push_back(i);
        }
    }
    return tmp;
}

template<typename T>
static std::vector<size_t> evalPos(const T& step, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    std::size_t idx{0};
    auto cmp = [&filter](const Vec& at){
        size_t dir = (filter.pos & filter.DIR_MASK) >> 3;
        if(filter.pos & filter.P_LT){
            return at[dir] < filter.posVal;
        }
        return at[dir] > filter.posVal;
    };
    for(const auto& at: step.asFmt(static_cast<AtomFmt>(filter.pos & filter.FMT_MASK))){
        if(cmp(at.coord)){
            tmp.push_back(idx);
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static std::vector<size_t> evalCoord(const T& step, const SelectionFilter& filter){
    std::vector<size_t> tmp;
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
            tmp.push_back(i);
        }
    }
    return tmp;
}

template<typename T>
static std::vector<size_t> invertSel(const T& step, const std::vector<size_t>& in){
    std::vector<size_t> out;
    const auto nat = step.getNat();
    out.reserve(nat);
    for(size_t i=0; i<nat; ++i){
        if(std::find(in.begin(), in.end(), i) == in.end()){
            out.push_back(i);
        }
    }
    return out;
}

template<typename T>
static std::vector<size_t> evalSubFilter(const T& step,
                                  const SelectionFilter& filter,
                                  SelectionFilter& subfilter,
                                  const std::vector<size_t>& parent){
    using Op = SelectionFilter::Op;
    std::vector<size_t> child = evalFilter(step, subfilter);
    std::vector<size_t> tmp(parent.size()+child.size());
    std::vector<size_t>::iterator it;
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
std::vector<size_t> evalFilter(const T& step, SelectionFilter& filter)
{
    std::vector<size_t> tmp;
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

/*
 * Selection container
 *
 * contains indices of selected atoms in Step-like
 */
template<typename T>
struct AtomSelection{
    std::vector<size_t> indices;
    T*                  step;
    SelectionFilter     filter;

    template<typename A>
    class AtomSelIterator: private T::iterator
    {
    private:
        using Base = typename T::iterator;
    public:
        AtomSelIterator(std::shared_ptr<AtomSelection> selection,
                        AtomFmt, size_t idx)
        //TODO: introduce a terminal-object (when c++17 is ready?)
            : Base{selection->step->begin()},
              selection{selection}, idx{idx}
        {
            size_t diff = selection->indices.empty() ? 0 : selection->indices[idx];
            Base::operator+=(diff);
        }
        AtomSelIterator& operator++(){
            return operator+=(1);
        }
        AtomSelIterator& operator+=(size_t i){
            idx += i;
            auto diff = selection->indices[idx] - selection->indices[idx-i];
            Base::operator+=(diff);
            return *this;
        }
        AtomSelIterator operator+(size_t i){
            AtomSelIterator copy{*this};
            return copy+=i;
        }
        A&  operator*() const {
            return Base::operator*();
        }
        A*  operator->() const {
            return Base::operator->();
        }
        bool    operator==(const AtomSelIterator& rhs) const noexcept{
            return (selection == rhs.selection) && (idx == rhs.idx);
        }
        bool    operator!=(const AtomSelIterator& rhs) const noexcept{
            return !(*this == rhs);
        }
        size_t getIdx() const noexcept{
            return idx;
        }
    private:
        std::shared_ptr<AtomSelection> selection;
        size_t idx;
    };

    using iterator = AtomSelIterator<Atom>;
    using constIterator = AtomSelIterator<const Atom>;
};

/*
 * Basic Selection-class template
 *
 * Instantiation of Bond- and Cell-getters with AtomSelection as Atom-source
 * Shall be instantiated with StepConst or StepMutable depending on desired const-ness,
 * and Step-like with desired target AtomSource
 */

template<template<typename> class B, typename T>
class SelectionBase: public B<AtomSelection<B<T>>>
{
private:
    using SelStep = B<T>;
    using Source = AtomSelection<SelStep>;
    using Base = B<Source>;
public:
    SelectionBase(std::shared_ptr<PseMap> p, AtomFmt f,
                  const SelStep* s, SelectionFilter sf,
                  std::shared_ptr<BondList> b,
                  std::shared_ptr<CellData> c, std::shared_ptr<std::string> co)
        : Base{p, f, std::make_shared<Source>(), b, c, co}
    {
        this->atoms->step = const_cast<SelStep*>(s);
        setFilter(sf);
        this->evaluateCache();
    }
    SelectionBase(std::shared_ptr<PseMap> p, AtomFmt f,
                  const SelStep* s, std::string sf,
                  std::shared_ptr<BondList> b,
                  std::shared_ptr<CellData> c, std::shared_ptr<std::string> co)
        : Base{p, f, std::make_shared<Source>(), b, c, co}
    {
        this->atoms->step = const_cast<SelStep*>(s);
        setFilter(sf);
        this->evaluateCache();
    }
    SelectionBase(std::shared_ptr<PseMap> p, AtomFmt f,
                  std::shared_ptr<Source> a,
                  std::shared_ptr<BondList> b,
                  std::shared_ptr<CellData> c, std::shared_ptr<std::string> co)
        : Base{p, f, a, b, c, co}
    {}
    SelectionBase(const SelectionBase& s)
        : Base{s.pse, s.at_fmt,
               std::make_shared<Source>(*s.atoms),
               std::make_shared<BondList>(*s.bonds),
               s.cell,
               std::make_shared<std::string>(*s.comment)}
    {}
    SelectionBase& operator=(const SelectionBase& s)
    {
        this->pse = s.pse;
        this->at_fmt = s.at_fmt;
        *this->bonds = *s.bonds;
        this->cell = s.cell;
        *this->comment = *s.comment;
        *this->atoms = *s.atoms;
        this->evaluateCache();
        return *this;
    }

    SelectionBase asFmt(AtomFmt tgt) // hides StepMutable::asFmt
    {
        auto tmp = SelectionBase{this->pse, tgt, this->atoms,
                this->bonds, this->cell, this->comment};
        tmp.evaluateCache();
        return tmp;
    }
    using StepConst<Source>::asFmt;

    const std::vector<size_t>& getIndices() const noexcept
    {
        return this->atoms->indices;
    }

    const SelectionFilter& getFilter() const noexcept
    {
        return this->atoms->filter;
    }
    void setFilter(std::string filter)
    {
        auto fs = std::stringstream{filter};
        fs >> this->atoms->filter;
    }
    void setFilter(SelectionFilter filter)
    {
        this->atoms->filter = std::move(filter);
    }
    void evaluateFilter() const
    {
        this->atoms->indices = evalFilter(this->atoms->step, this->atoms->filter);
    }
};

}

#endif // LIBVIPSTER_STEPSEL_H
