#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include "stepbase.h"
#include "atomcontainers.h"

namespace Vipster {
// Forward declarations
class Step;
template<typename T>
class SelectionProper;
using StepSelection = SelectionProper<Step>;
using StepSelConst = SelectionProper<const Step>;
template<typename T>
class SelectionFormatter;
using StepSelFormatter = SelectionFormatter<Step>;

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
std::vector<size_t> evalFilter(const Step& step, SelectionFilter& filter);

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
 * Basic Selection-class template
 *
 * Instantiation of Bond- and Cell-getters with AtomSelection as Atom-source
 * Shall be instanced with `Step` or `const Step` as template argument
 */
template<typename T>
class SelectionBase: public StepBase<SelectionBase<T>>
{
public:
    SelectionBase& operator=(const SelectionBase& s)
    {
        this->pse = s.pse;
        this->at_fmt = s.at_fmt;
        *this->bonds = *s.bonds;
        this->cell = s.cell;
        *this->comment = *s.comment;
        *this->selection = *s.selection;
        *this->filter = *s.filter;
        this->step = s.step;
        return *this;
    }
    // TODO: missing constructors?
    // TODO: direct Prop-getters?
    const std::vector<size_t>& getIndices() const noexcept
    {
        evaluateCache();
        return selection->indices;
    }

    void        evaluateCache() const override
    {
        // make sure that caches are clean, for pos even the needed formatted cache
        auto fmt = (filter->mode == SelectionFilter::Mode::Pos) ?
                    static_cast<AtomFmt>(filter->pos & SelectionFilter::FMT_MASK) :
                    this->at_fmt;
        step.asFmt(fmt).evaluateCache();
        if(filter->op & SelectionFilter::UPDATE){
            selection->indices = evalFilter(step, *filter);
        }
    }
    Vec   getCenter(CdmFmt fmt, bool com=false) const noexcept override
    {
        return step.getCenter(fmt, com);
    }
    const SelectionFilter& getFilter() const noexcept
    {
        return *filter;
    }
    void setFilter(std::string filter)
    {
        auto fs = std::stringstream{filter};
        fs >> *this->filter;
    }
    void setFilter(SelectionFilter filter)
    {
        *this->filter = std::move(filter);
    }
    void evaluateFilter() const
    {
        selection->indices = evalFilter(step, *filter);
    }

    // Atoms
    size_t          getNat() const noexcept
    {
        evaluateCache();
        return selection->indices.size();
    }
    using           iterator = AtomSelIterator<Atom>;
    using           constIterator = AtomSelIterator<const Atom>;
    Atom            operator[](size_t i)
    {
        evaluateCache();
        return *iterator{selection, this->at_fmt, i};
    }
    iterator        begin() noexcept
    {
        evaluateCache();
        return iterator{selection, this->at_fmt, 0};
    }
    constIterator   begin() const noexcept
    {
        evaluateCache();
        return constIterator{selection, this->at_fmt, 0};
    }
    constIterator   cbegin() const noexcept
    {
        evaluateCache();
        return constIterator{selection, this->at_fmt, 0};
    }
    iterator        end() noexcept
    {
        evaluateCache();
        return iterator{selection, this->at_fmt, selection->indices.size()};
    }
    constIterator   end() const noexcept
    {
        evaluateCache();
        return constIterator{selection, this->at_fmt, selection->indices.size()};
    }
    constIterator   cend() const noexcept
    {
        evaluateCache();
        return constIterator{selection, this->at_fmt, selection->indices.size()};
    }
    void delAtoms()
    {
        auto& idx = selection->indices;
        for(auto it=idx.rbegin(); it!= idx.rend(); ++it){
            step.delAtom(*it);
        }
        setFilter(SelectionFilter{});
    }
protected:
    SelectionBase(T& step)
        : StepBase<SelectionBase<T>>{step.pse,
                   step.at_fmt,
                   std::make_shared<BondList>(),
                   step.cell,
                   std::make_shared<std::string>(*step.comment)},
          selection{std::make_shared<AtomSelection>(
                        AtomSelection{std::vector<size_t>{}, step.atoms})},
          filter{std::make_shared<SelectionFilter>()},
          step{step}
    {}
    SelectionBase(const SelectionBase& s)
        : StepBase<SelectionBase<T>>{s.pse,
                   s.at_fmt,
                   std::make_shared<BondList>(*s.bonds),
                   s.cell,
                   std::make_shared<std::string>(*s.comment)},
          selection{std::make_shared<AtomSelection>(*s.selection)},
          filter{std::make_shared<SelectionFilter>(*s.filter)},
          step{s.step}
    {}
    SelectionBase(std::shared_ptr<PseMap> p, AtomFmt f, std::shared_ptr<BondList> b,
         std::shared_ptr<CellData> c, std::shared_ptr<std::string> s,
         std::shared_ptr<AtomSelection> a, std::shared_ptr<SelectionFilter> sf,
         T& step)
        : StepBase<SelectionBase<T>>{std::move(p), std::move(f),
                   std::move(b), std::move(c), std::move(s)},
          selection{std::move(a)}, filter{std::move(sf)}, step{step}
    {}
    std::shared_ptr<AtomSelection>      selection;
    std::shared_ptr<SelectionFilter>    filter;
    T&                                  step;
};

template<typename T>
class SelectionFormatter: public SelectionBase<T>
{
public:
    SelectionFormatter(SelectionProper<T>& sel, AtomFmt at_fmt)
        : SelectionBase<T>{sel.pse,
                           at_fmt,
                           sel.bonds,
                           sel.cell,
                           sel.comment,
                           sel.selection,
                           sel.filter,
                           sel.step},
          sel{sel}
    {}
    SelectionBase<T>&       asFmt(AtomFmt fmt) override
    {
        return sel.asFmt(fmt);
    }
    const SelectionBase<T>& asFmt(AtomFmt fmt) const override
    {
        return sel.asFmt(fmt);
    }
private:
    SelectionProper<T>& sel;
};

template<typename T>
class SelectionProper: public SelectionBase<T>
{
    friend class SelectionFormatter<T>;
public:
    SelectionProper(T& step)
        : SelectionBase<T>{step}
    {}
    SelectionProper(T& step, std::string filter)
        : SelectionProper{step}
    {
        this->setFilter(filter);
    }
    SelectionProper(T& step, SelectionFilter filter)
        : SelectionProper{step}
    {
        this->setFilter(filter);
    }
    SelectionProper(const SelectionProper& s)
        :SelectionBase<T>{s}
    {}
    SelectionProper& operator=(const SelectionProper& s){
        *static_cast<SelectionBase<T>*>(this) = s;
        return *this;
    }
    void    setFmt(AtomFmt at_fmt)
    {
        this->at_fmt = at_fmt;
    }
    SelectionFormatter<T>   asBohr{*this, AtomFmt::Bohr};
    SelectionFormatter<T>   asAngstrom{*this, AtomFmt::Angstrom};
    SelectionFormatter<T>   asCrystal{*this, AtomFmt::Crystal};
    SelectionFormatter<T>   asAlat{*this, AtomFmt::Alat};
    SelectionBase<T>&       asFmt(AtomFmt fmt) override
    {
        switch(fmt){
        case AtomFmt::Bohr:
            return asBohr;
        case AtomFmt::Angstrom:
            return asAngstrom;
        case AtomFmt::Crystal:
            return asCrystal;
        default:
            return asAlat;
        }
    }
    const SelectionBase<T>& asFmt(AtomFmt fmt) const override
    {
        return const_cast<SelectionProper<T>*>(this)->asFmt(fmt);
    }
};

}

#endif // LIBVIPSTER_STEPSEL_H
