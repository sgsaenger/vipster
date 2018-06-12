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
struct SelectionFilter{
    SelectionFilter() = default;
    SelectionFilter(const SelectionFilter& f) {
        *this = f;
    }
    SelectionFilter& operator=(const SelectionFilter& f){
        mode = f.mode; op = f.op;
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
    enum class Mode:uint8_t{Index, Type, Coord, Pos, Group};
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
    uint8_t op;
    uint8_t pos;
    uint8_t coord;
    float posVal;
    size_t coordVal;
    std::set<size_t> indices;
    std::set<std::string> types;
    std::unique_ptr<SelectionFilter> groupfilter{nullptr};
    std::unique_ptr<SelectionFilter> subfilter{nullptr};
};

std::ostream& operator<<(std::ostream& os, const SelectionFilter& filter);
std::istream& operator>>(std::istream& is, SelectionFilter& filter);
std::vector<size_t> evalFilter(const Step& step, SelectionFilter& filter);

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
        *this->cell = *s.cell;
        *this->comment = *s.comment;
        *this->selection = *s.selection;
        *this->filter = *s.filter;
        this->step = s.step;
        return *this;
    }
    // TODO: missing constructors?
    // TODO: direct Prop-getters?

    void        evaluateCache() const override
    {
        using Mode = SelectionFilter::Mode;
        bool updateNeeded{static_cast<bool>(filter->op & SelectionFilter::UPDATE)};
        // Type or Index need to be reset when types have been changed
        if(!updateNeeded &&
                ((filter->mode == Mode::Type) ||
                 (filter->mode == Mode::Index))){
            if(step.atoms->name_changed)
                updateNeeded = true;
        }
        // Pos or Index need to be reset when coords have been changed
        if(!updateNeeded &&
                ((filter->mode == Mode::Pos) ||
                 (filter->mode == Mode::Index))){
            if(std::any_of(step.atoms->coord_changed.begin(), // any coordinates changed,
                           step.atoms->coord_changed.end(),   // no caches evaluated
                           [](bool b){return b;}) ||
                step.atoms->coord_outdated[filter->pos &      // something changed,
                     SelectionFilter::FMT_MASK]){             // relevant cache still outdated
                updateNeeded = true;
            }
        }
        // make sure that caches are clean, for pos even the needed formatted cache
        auto fmt = (filter->mode == Mode::Pos) ?
                    static_cast<AtomFmt>(filter->pos & SelectionFilter::FMT_MASK) :
                    this->at_fmt;
        step.asFmt(fmt).evaluateCache();
        // Coord needs to be reset when bonds (possibly) change
        if(!updateNeeded && (filter->mode == Mode::Coord)){
            if(step.bonds->outdated){
                updateNeeded = true;
            }
        }
        if(updateNeeded){
            selection->indices = evalFilter(step, *filter);
        }
    }
    std::string getFilter() const noexcept
    {
        std::stringstream ss{};
        ss << *filter;
        return ss.str();
    }
    void setFilter(std::string filter)
    {
        auto fs = std::stringstream{filter};
        fs >> *this->filter;
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
protected:
    SelectionBase(T& step)
        : StepBase<SelectionBase<T>>{step.pse,
                   step.at_fmt,
                   std::make_shared<BondList>(),
                   std::make_shared<CellData>(*step.cell),
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
                   std::make_shared<CellData>(*s.cell),
                   std::make_shared<std::string>(*s.comment)},
          selection{std::make_shared<AtomSelection>(*s.selection)},
          filter{std::make_shared<SelectionFilter>(*s.filter)},
          step{s.step}
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
        : SelectionBase<T>{sel.step}, sel{sel}
    {
        this->at_fmt = at_fmt;
    }
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
        case AtomFmt::Alat:
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
