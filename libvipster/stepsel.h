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

using FilterPair = std::pair<size_t, std::vector<SizeVec>>;

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
             X=0x0, Y=0x4, Z=0x8, DIR_MASK=0xC, // 2  bits for space direction
             P_GT=0x0, P_LT=0x10, P_CMP_MASK=0x10,  // 1 bit for comp direction
            };
    enum Coord{C_GT=0x0, C_EQ=0x1, C_LT=0x2, C_CMP_MASK=0x3};
    Mode mode;
    uint8_t op{Op::UPDATE};
    uint8_t pos;
    uint8_t coord;
    float posVal;
    size_t coordVal;
    std::map<size_t, std::vector<SizeVec>> indices;
    std::set<std::string> types;
    std::unique_ptr<SelectionFilter> groupfilter{nullptr};
    std::unique_ptr<SelectionFilter> subfilter{nullptr};
};

/*
 * recursively evaluate the filters
 */
template<typename T>
static std::vector<FilterPair> evalType(const T& step, const SelectionFilter& filter){
    std::vector<FilterPair> tmp;
    size_t idx{0};
    for(const auto& at: step){
        for(const auto& type: filter.types){
            if(at.name == type){
                tmp.emplace_back(idx, std::vector{SizeVec{}});
                break;
            }
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static std::vector<FilterPair> evalIdx(const T& step, const SelectionFilter& filter){
    std::vector<FilterPair> tmp;
    auto nat = step.getNat();
    for(const auto& p:filter.indices){
        if(p.first<nat){
            tmp.push_back(p);
        }
    }
    return tmp;
}

template<typename T>
static std::vector<FilterPair> evalPos(const T& step, const SelectionFilter& filter){
    std::vector<FilterPair> tmp;
    std::size_t idx{0};
    auto cmp = [&filter](const Vec& at){
        size_t dir = (filter.pos & filter.DIR_MASK) >> 2;
        if(filter.pos & filter.P_LT){
            return at[dir] < filter.posVal;
        }
        return at[dir] > filter.posVal;
    };
    for(const auto& at: step.asFmt(static_cast<AtomFmt>(filter.pos & filter.FMT_MASK))){
        if(cmp(at.coord)){
            tmp.emplace_back(idx, std::vector{SizeVec{}});
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static std::vector<FilterPair> evalCoord(const T& step, const SelectionFilter& filter){
    std::vector<FilterPair> tmp;
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
            tmp.emplace_back(i, std::vector{SizeVec{}});
        }
    }
    return tmp;
}

template<typename T>
static std::vector<FilterPair> invertSel(const T& step, const std::vector<FilterPair>& in){
    std::vector<FilterPair> out;
    const auto nat = step.getNat();
    out.reserve(nat);
    for(size_t i=0; i<nat; ++i){
        if(std::find_if(in.begin(), in.end(),
                        [&i](const auto& pair){
                            return pair.first == i;
                        }) == in.end()){
            out.emplace_back(i, std::vector{SizeVec{}});
        }
    }
    return out;
}

template<typename T>
static std::vector<FilterPair> evalSubFilter(const T& step,
                                  const SelectionFilter& filter,
                                  SelectionFilter& subfilter,
                                  const std::vector<FilterPair>& parent){
    using Op = SelectionFilter::Op;
    std::vector<FilterPair> child = evalFilter(step, subfilter);
    std::vector<FilterPair> tmp(parent.size()+child.size());
    std::vector<FilterPair>::iterator it;
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
std::vector<FilterPair> evalFilter(const T& step, SelectionFilter& filter)
{
    std::vector<FilterPair> tmp;
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
    std::vector<FilterPair> indices;
    // TODO: molecule is not performance-critical. use shared-ptrs here and in molecule?
    T*                  step;
    SelectionFilter     filter;

    size_t getNat() const noexcept
    {
        return indices.size();
    }

    void evaluateCache(const StepConst<AtomSelection<T>>& step)
    {
        auto fmt = (filter.mode == SelectionFilter::Mode::Pos) ?
                    static_cast<AtomFmt>(filter.pos & SelectionFilter::FMT_MASK) :
                    step.at_fmt;
        this->step->asFmt(fmt).evaluateCache();
        if(filter.op & SelectionFilter::UPDATE){
            indices = evalFilter(*this->step, filter);
        }
    }

    template<typename B>
    class AtomSelIterator;
    using iterator = AtomSelIterator<typename T::iterator>;
    using const_iterator = AtomSelIterator<typename T::const_iterator>;

    template<typename B>
    class AtomSelIterator: private B
    {
    public:
        using difference_type = ptrdiff_t;
        using value_type = typename B::value_type;
        using reference = typename B::reference;
        using pointer = typename B::pointer;
        using iterator_category = std::bidirectional_iterator_tag;
        AtomSelIterator(std::shared_ptr<AtomSelection> selection,
                        AtomFmt fmt, size_t idx)
        //TODO: introduce a terminal-object (when c++17 is ready?)
            : B{selection->step->asFmt(fmt).begin()},
              selection{selection}, idx{idx}
        {
            size_t diff = selection->indices.empty() ? 0 : selection->indices[idx].first;
            B::operator+=(diff);
        }
        // allow iterator to const_iterator conversion
        template <typename U, typename R=B,
                  typename = typename std::enable_if<std::is_same<typename T::const_iterator, R>::value>::type>
        AtomSelIterator(const AtomSelIterator<U>& it)
            : B{*it}, selection{it.selection}, idx{it.idx}
        {}
        AtomSelIterator& operator++(){
            return operator+=(1);
        }
        AtomSelIterator& operator--(){
            return operator+=(-1);
        }
        AtomSelIterator& operator+=(long i){
            idx += i;
            auto diff = selection->indices[idx].first - selection->indices[idx-i].first;
            B::operator+=(diff);
            return *this;
        }
        AtomSelIterator operator+(long i){
            AtomSelIterator copy{*this};
            return copy+=i;
        }
        AtomSelIterator operator-(long i){
            AtomSelIterator copy{*this};
            return copy+=i;
        }
        reference  operator*() const {
            return B::operator*();
        }
        pointer  operator->() const {
            return B::operator->();
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
        const FilterPair& getFilterPair() const{
            return selection->indices[idx];
        }
    private:
        std::shared_ptr<AtomSelection> selection;
        size_t idx;
    };
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
    template<template<typename> class Bt, typename Tt> friend class SelectionBase;
public:
    using SelSource = T;
    using SelStep = B<SelSource>;
    using Source = AtomSelection<SelStep>;
    using Base = B<Source>;
    SelectionBase(std::shared_ptr<PseMap> p, AtomFmt f,
                  const SelStep* s, SelectionFilter sf,
                  std::shared_ptr<CellData> c,
                  std::shared_ptr<std::string> co)
        : Base{p, f, std::make_shared<Source>(), std::make_shared<BondList>(), c, co}
    {
        this->atoms->step = const_cast<SelStep*>(s);
        setFilter(sf);
        this->evaluateCache();
    }
    SelectionBase(std::shared_ptr<PseMap> p, AtomFmt f,
                  const SelStep* s, std::string sf,
                  std::shared_ptr<CellData> c,
                  std::shared_ptr<std::string> co)
        : Base{p, f, std::make_shared<Source>(), std::make_shared<BondList>(), c, co}
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
               s.cell, s.comment}
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
    // mutable to const selection
    template <template<typename> class U, typename V,
              typename R=SelSource,
              typename = typename std::enable_if<std::is_same<StepConst<T>, R>::value>::type, // only allow for const_selection
              typename = typename std::enable_if<std::is_same<V, T>::value>::type> // if source is of same type
    SelectionBase(const SelectionBase<U,V>& sel)
        : Base{sel.pse, sel.at_fmt,
               std::make_shared<Source>(*sel.atoms),
               std::make_shared<BondList>(),
               sel.cell, sel.comment}
    {}
    // remove one level of selection-indirection
    template <template<typename> class U, class V,
              typename = typename std::enable_if<std::is_same<V, SelStep>::value>::type>
    SelectionBase(const SelectionBase<U,AtomSelection<V>>& sel)
        : Base{sel.pse, sel.at_fmt,
               std::make_shared<Source>(),
               std::make_shared<BondList>(),
               sel.cell, sel.comment}
    {
        auto& tmp = sel.getAtoms();
        auto& indices = tmp.indices;
        auto& remote_source = tmp.step->getAtoms();
        this->atoms->step = remote_source.step;
        SelectionFilter fil{};
        fil.mode = SelectionFilter::Mode::Index;
        for(auto& i: indices){
            fil.indices.insert(remote_source.indices[i.first]);
        }
        this->evaluateCache();
    }

    SelectionBase asFmt(AtomFmt tgt) // hides StepMutable::asFmt
    {
        auto tmp = SelectionBase{this->pse, tgt, this->atoms,
                this->bonds, this->cell, this->comment};
        tmp.evaluateCache();
        return tmp;
    }
    using StepConst<Source>::asFmt;

    const std::vector<FilterPair>& getIndices() const noexcept
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
        "<li><b><tt>pos</tt></b>: relative position (x,y,z to choose axis, a,b,c,d to choose format (Ångstrom, Bohr, Crystal, Alat, respectively))"
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

#endif // LIBVIPSTER_STEPSEL_H
