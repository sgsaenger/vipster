#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include "step.h"

namespace Vipster {
/*
 * Wrap multiple filter criteria without polymorphism
 * also, represents differently coupled filter chains
 * Target grammar:
 * Keywords are case insensitive, whereas atom-types are case-sensitive
 *
 * Filter ::= Criterion, [("(", Coupling, Filter, ")")];
 * Criterion ::= ["not "], (TypeCrit | IdxCrit | PosCrit | CoordCrit);
 * Coupling ::= ["!"], ("|" | "&" | "^");
 *
 * TypeCrit ::= "type ", (Type | TypeList);
 * TypeList ::= "[", Type, {(" ", Type)}, "]";
 * Type ::= NonWhiteSpace, {NonWhiteSpace};
 *
 * IdxCrit ::= "index " , IdxList;
 * IdxList ::= ( "[", IdxRange, {(" " IdxRange)}, "]") | IdxRange;
 * IdxRange ::= ( Integer, "-", Integer) | Integer;
 *
 * PosCrit ::= "pos ", Direction, CompOp, Float;
 * Direction ::= "x" | "y" | "z";
 * CompOp ::= ">" | "<";
 *
 * CoordCrit ::= "coord ", CompEqOp, Integer;
 * CompEqOp ::= "=" | CompOp;
 */
struct SelectionFilter{
    SelectionFilter() = default;
    SelectionFilter(const SelectionFilter& f)
        :mode{f.mode}, op{f.op},
          indices{f.indices},
          types{f.types},
          subfilter{std::make_unique<SelectionFilter>(*f.subfilter)}
    {}
    enum class Mode{Index, Type, Coord, Pos};
    enum Op{NONE=0x0, NOT=0x1, // first bit negates own op
            PAIR=0x2, NOT_PAIR=0x4, // second bit activates coupling, third bit negates
            AND=0x2, NAND=0x6,
            OR=0xA, NOR=0xE,
            XOR=0x12, XNOR=0x16,
            UPDATE=0x80};
    Mode mode;
    unsigned int op;
    std::vector<size_t> indices;
    std::vector<std::string> types;
    std::unique_ptr<SelectionFilter> subfilter{nullptr};
};

std::stringstream& operator<<(std::stringstream& os, const SelectionFilter& filter)
{
    using Op = SelectionFilter::Op;
    using Mode = SelectionFilter::Mode;
    if(filter.op & Op::NOT){
        os << "not ";
    }
    switch(filter.mode){
    case Mode::Index:
        os << "index ";
        // TODO
        break;
    case Mode::Coord:
        os << "coord ";
        // TODO
        break;
    case Mode::Pos:
        os << "pos ";
        // TODO
        break;
    case Mode::Type:
        os << "type ";
        if(filter.types.size() == 1){
            os << filter.types.back();
        }else{
            os << "[ ";
            for(const auto& t:filter.types){
                os << t << " ";
            }
            os << ']';
        }
        break;
    }
    if(filter.op & Op::PAIR){
        os << '(';
        if(filter.op & Op::NOT_PAIR){
            os << '!';
        }
        if(filter.op & Op::OR){
            os << '|';
        }else if(filter.op & Op::XOR){
            os << '^';
        }else{
            os << '&';
        }
        os << *filter.subfilter;
        os << ')';
    }
    return os;
}

std::stringstream& operator>>(std::stringstream& is, SelectionFilter& filter){
    using Op = SelectionFilter::Op;
    using Mode = SelectionFilter::Mode;
    std::string token;
    filter.op = 0;
    is >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if(token == "not"){
        filter.op |= Op::NOT;
        is >> token;
        std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    }
    if(token == "type"){
        filter.mode = Mode::Type;
        is >> token;
        if(token[0] == '['){
            if(token.size() > 1){
                filter.types.push_back(token.substr(1));
            }
            while(is.good()){
                is >> token;
                if(token.back() == ']'){
                    token.pop_back();
                    if(!token.empty()){
                        filter.types.push_back(token);
                    }
                    break;
                }else{
                    filter.types.push_back(token);
                }
            }
        }else{
            is >> token;
            filter.types.push_back(token);
        }
    }else if(token == "coord"){
        filter.mode = Mode::Coord;
        // TODO
    }else if(token == "pos"){
        filter.mode = Mode::Pos;
        // TODO
    }else if(token == "index"){
        filter.mode = Mode::Index;
        // TODO
    }
    auto ctoken = is.peek();
    if(ctoken == '('){
        is.get();
        ctoken = is.get();
        if(ctoken == '!'){
            filter.op |= Op::NOT_PAIR;
            ctoken = is.get();
        }
        if(ctoken == '|'){
            filter.op |= Op::OR;
        }else if(ctoken == '&'){
            filter.op |= Op::AND;
        }else if(ctoken == '^'){
            filter.op |= Op::XOR;
        }else{
            throw Error("Unknown coupling operator "+std::string{1, static_cast<char>(ctoken)});
        }
        filter.subfilter = std::make_unique<SelectionFilter>();
        is >> *filter.subfilter;
        ctoken = is.get();
        if(ctoken != ')'){
            throw Error("Unterminated nested filter");
        }
    }
    filter.op |= Op::UPDATE;
    return is;
}

/*
 * Selection container
 *
 * contains indices of selected atoms in AtomList
 */
struct AtomSelection{
    std::vector<size_t> indices;
    std::shared_ptr<AtomList> atoms;
};

/*
 * recursively evaluate the filters
 */
std::vector<size_t> evalFilter(const AtomList& atoms, SelectionFilter& filter)
{
    std::vector<size_t> tmp;
    if(filter.mode == SelectionFilter::Mode::Type){
        size_t idx{0};
        for(const auto& at: atoms.names){
            for(const auto& type: filter.types){
                if(at == type){
                    tmp.push_back(idx);
                    break;
                }
            }
            ++idx;
        }
    }
    return tmp;
}

/*
 * Iterator for Atom selection
 *
 * dereferences selection-indices
 */
template<typename T>
class AtomSelIterator: private T
{
public:
    AtomSelIterator(const std::shared_ptr<AtomSelection> &selection,
                    AtomFmt fmt, size_t idx)
        : T{&selection->atoms->coordinates[static_cast<size_t>(fmt)][selection->indices[idx]],
            &selection->atoms->coord_changed[static_cast<size_t>(fmt)],
            &selection->atoms->names[selection->indices[idx]],
            &selection->atoms->name_changed,
            &selection->atoms->charges[selection->indices[idx]],
            &selection->atoms->properties[selection->indices[idx]],
            &selection->atoms->pse[selection->indices[idx]],
            &selection->atoms->prop_changed},
          selection{selection}, fmt{fmt}, idx{idx}
    {}
    AtomSelIterator& operator++(){
        ++idx;
        auto diff = selection->indices[idx] - selection->indices[idx-1];
        this->coord_ptr += diff;
        this->name_ptr += diff;
        this->charge_ptr += diff;
        this->prop_ptr += diff;
        this->pse_ptr += diff;
        return *this;
    }
    AtomSelIterator& operator+=(size_t i){
        idx += i;
        auto diff = selection->indices[idx] - selection->indices[idx-i];
        this->coord_ptr += diff;
        this->name_ptr += diff;
        this->charge_ptr += diff;
        this->prop_ptr += diff;
        this->pse_ptr += diff;
        return *this;
    }
    AtomSelIterator operator+(size_t i){
        AtomSelIterator copy{*this};
        return copy+=i;
    }
    T&  operator*() const {
        return static_cast<T&>(*const_cast<AtomSelIterator*>(this));
    }
    T*  operator->() const {
        return &(operator*());
    }
    bool    operator==(const AtomSelIterator& rhs) const noexcept{
        return (selection == rhs.selection) && (fmt == rhs.fmt) && (idx == rhs.idx);
    }
    bool    operator!=(const AtomSelIterator& rhs) const noexcept{
        return !(*this == rhs);
    }
    size_t getIdx() const noexcept{
        return idx;
    }
private:
    std::shared_ptr<AtomSelection> selection;
    AtomFmt fmt;
    size_t idx;
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
        pse = s.pse;
        this->at_fmt = s.at_fmt;
        this->atoms = std::make_shared<AtomList>(*s.atoms);
        this->bonds = std::make_shared<BondList>(*s.bonds);
        this->cell = std::make_shared<CellData>(*s.cell);
        this->comment = std::make_shared<std::string>(*s.comment);
        this->selection = std::make_shared<AtomSelection>(*s.selection);
        this->filter = std::make_shared<SelectionFilter>(*s.filter);
        this->step = s.step;
        return *this;
    }
    //TODO: missing constructors?
    // TODO: direct Prop-getters?

    void        evaluateCache() const override
    {
        bool updateNeeded{static_cast<bool>(filter->op & SelectionFilter::UPDATE)};
        if(!updateNeeded){
            switch(filter->mode){
            case SelectionFilter::Mode::Type:
                if(step.atoms->name_changed)
                    updateNeeded = true;
                break;
            case SelectionFilter::Mode::Coord:
                // TODO
                break;
            case SelectionFilter::Mode::Index:
                // TODO
                break;
            case SelectionFilter::Mode::Pos:
                // TODO
                break;
            }
        }
        step.evaluateCache();
        if(updateNeeded){
            selection->indices = evalFilter(*selection->atoms, *filter);
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
