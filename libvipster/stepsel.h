#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include "step.h"

namespace Vipster {

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
 * Wrap multiple filter criteria without polymorphism
 * also, represents differently coupled filter chains
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
    SelectionFilter(const SelectionFilter& f)
        :mode{f.mode}, op{f.op},
          indices{f.indices},
          types{f.types},
          subfilter{std::make_unique<SelectionFilter>(*f.subfilter)}
    {}
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

void writeType(std::ostream& os, const SelectionFilter& filter){
    os << "type ";
    if(filter.types.size() == 1){
        os << *filter.types.cbegin();
    }else{
        os << "[ ";
        for(const auto& t:filter.types){
            os << t << " ";
        }
        os << ']';
    }
}

void writeCoord(std::ostream& os, const SelectionFilter& filter){
    using Coord = SelectionFilter::Coord;
    os << "coord ";
    switch(filter.coord & Coord::C_CMP_MASK){
    case Coord::C_GT:
        os << '>';
        break;
    case Coord::C_EQ:
        os << '=';
        break;
    case Coord::C_LT:
        os << '<';
        break;
    }
    os << filter.coordVal;
}

void writePos(std::ostream& os, const SelectionFilter& filter){
    using Pos = SelectionFilter::Pos;
    os << "pos ";
    const std::map<uint8_t, char> posToC{
        {Pos::X, 'x'}, {Pos::Y, 'y'}, {Pos::Z, 'z'},
        {Pos::ANG, 'a'}, {Pos::BOHR, 'b'},
        {Pos::CRYS, 'c'}, {Pos::CDM, 'd'},
        {Pos::P_LT, '<'}, {Pos::P_GT, '>'}
    };
    os << posToC.at(filter.pos & Pos::DIR_MASK);
    os << posToC.at(filter.pos & Pos::FMT_MASK);
    os << posToC.at(filter.pos & Pos::P_CMP_MASK);
    os << filter.posVal;
}

void writeIdx(std::ostream& os, const SelectionFilter& filter){
    os << "index ";
    if(filter.indices.size() == 1){
        os << *filter.indices.cbegin();
    }else{
        os << "[ ";
        for(const auto& i:filter.indices){
            os << i << " ";
        }
        os << ']';
    }
}

std::ostream& operator<<(std::ostream& os, const SelectionFilter& filter)
{
    using Op = SelectionFilter::Op;
    using Mode = SelectionFilter::Mode;
    if(filter.op & Op::NOT){
        os << "not ";
    }
    switch(filter.mode){
    case Mode::Index:
        writeIdx(os, filter);
        break;
    case Mode::Coord:
        writeCoord(os, filter);
        break;
    case Mode::Pos:
        writePos(os, filter);
        break;
    case Mode::Type:
        writeType(os, filter);
        break;
    case Mode::Group:
        os << "( ";
        os << *filter.groupfilter;
        os << " )";
        break;
    }
    if(filter.op & Op::PAIR){
        os << ' ';
        if(filter.op & Op::NOT_PAIR){
            os << '!';
        }
        if((filter.op & Op::OR) == Op::OR){
            os << '|';
        }else if((filter.op & Op::XOR) == Op::XOR){
            os << '^';
        }else{
            os << '&';
        }
        os << ' ';
        os << *filter.subfilter;
    }
    return os;
}

void parseType(std::istream& is, SelectionFilter& filter){
    std::string token;
    filter.mode = SelectionFilter::Mode::Type;
    is >> token;
    if(token.front() == '['){
        if(token.size() > 1){
            filter.types.insert(token.substr(1));
        }
        while(is.good()){
            is >> token;
            if(token.back() == ')'){
                if(token.size() > 1){
                    token.pop_back();
                    is.putback(')');
                }else{
                    throw Error("Unterminated type list");
                }
            }
            if(token.back() == ']'){
                if(token.size() > 1){
                    token.pop_back();
                    filter.types.insert(token);
                }
                return;
            }else{
                filter.types.insert(token);
            }
        }
    }else{
        if(token.back() == ')'){
            if(token.size() > 1){
                token.pop_back();
                is.putback(')');
            }else{
                throw Error("Missing type");
            }
        }
        filter.types.insert(token);
    }
}

void parseCoord(std::istream& is, SelectionFilter& filter){
    using Coord = SelectionFilter::Coord;
    filter.mode = SelectionFilter::Mode::Coord;
    filter.coord = 0;
    char ctoken;
    is >> ctoken;
    switch(ctoken){
    case '<':
        filter.coord |= Coord::C_LT;
        break;
    case '>':
        filter.coord |= Coord::C_GT;
        break;
    case '=':
        filter.coord |= Coord::C_EQ;
        break;
    default:
        throw Error(std::string{"Unknown comparison "}+ctoken);
    }
    is >> filter.coordVal;
}

void parsePos(std::istream& is, SelectionFilter& filter){
    using Pos = SelectionFilter::Pos;
    filter.mode = SelectionFilter::Mode::Pos;
    filter.pos = 0;
    const std::map<char, Pos> cToPos{
        {'x', Pos::X}, {'y', Pos::Y}, {'z', Pos::Z},
        {'a', Pos::ANG}, {'b', Pos::BOHR},
        {'c', Pos::CRYS}, {'d', Pos::CDM},
        {'<', Pos::P_LT}, {'>', Pos::P_GT}
    };
    char ctoken;
    is >> ctoken;
    filter.pos |= cToPos.at(ctoken);
    is >> ctoken;
    filter.pos |= cToPos.at(ctoken);
    is >> ctoken;
    filter.pos |= cToPos.at(ctoken);
    is >> filter.posVal;
}

void parseIdx(std::istream& is, SelectionFilter& filter){
    filter.mode = SelectionFilter::Mode::Index;
    std::string token;
    is >> token;
    auto rangeParser = [&filter](const std::string &token){
        auto splitPos = token.find('-');
        if(splitPos != token.npos){
            auto left = std::stoul(token.substr(0, splitPos));
            auto right = std::stoul(token.substr(splitPos+1));
            for(auto i=std::min(left, right); i<=std::max(left,right);++i){
                filter.indices.insert(static_cast<size_t>(i));
            }
        }else{
            filter.indices.insert(std::stoul(token));
        }
    };
    if(token.front() == '['){
        if(token.size() > 1){
            rangeParser(token.substr(1));
        }
        while(is.good()){
            is >> token;
            if(token.back() == ')'){
                if(token.size() > 1){
                    token.pop_back();
                    is.putback(')');
                }else{
                    throw Error("Unterminated index list");
                }
            }
            if(token.back() == ']'){
                if(token.size() > 1){
                    token.pop_back();
                    rangeParser(token);
                }
                return;
            }else{
                rangeParser(token);
            }
        }
    }else{
        if(token.back() == ')'){
            if(token.size() > 1){
                token.pop_back();
                is.putback(')');
            }else{
                throw Error("Missing index");
            }
        }
        rangeParser(token);
    }
}

std::istream& operator>>(std::istream& is, SelectionFilter& filter){
    using Op = SelectionFilter::Op;
    std::string token;
    char c;
    filter.op = Op::UPDATE;
    is >> c;
    if(c == '('){
        // Parse subfilter (lower level, filter is parent)
        filter.mode = SelectionFilter::Mode::Group;
        filter.groupfilter = std::make_unique<SelectionFilter>();
        is >> *filter.groupfilter;
        is >> c;
        if(c != ')'){
            throw Error("Unterminated filter group");
        }
    }else{
        // Parse criterion
        is >> token;
        token.insert(token.begin(), c);
        std::transform(token.begin(), token.end(), token.begin(), ::tolower);
        if(token == "not"){
            filter.op |= Op::NOT;
            is >> token;
            std::transform(token.begin(), token.end(), token.begin(), ::tolower);
        }
        if(token == "type"){
            parseType(is, filter);
        }else if(token == "coord"){
            parseCoord(is, filter);
        }else if(token == "pos"){
            parsePos(is, filter);
        }else if(token == "index"){
            parseIdx(is, filter);
        }else{
            throw Error("Unknown selection operator: "+token);
        }
    }
    // Parse subfilter (same level)
    if(is.good() && (is >> token) && !token.empty()){
        if(token.front() == ')'){
            is.putback(')');
            return is;
        }
        if(token.size() > 2){
            throw Error("Unknown coupling operator "+token);
        }
        if(token.front() == '!'){
            filter.op |= Op::NOT_PAIR;
        }
        switch(token.back()){
        case '|':
            filter.op |= Op::OR;
            break;
        case '&':
            filter.op |= Op::AND;
            break;
        case '^':
            filter.op |= Op::XOR;
            break;
        default:
            throw Error("Unknown coupling operator "+token);
        }
        filter.subfilter = std::make_unique<SelectionFilter>();
        is >> *filter.subfilter;
    }
    return is;
}

/*
 * recursively evaluate the filters
 */
std::vector<size_t> evalType(const Step& step, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    size_t idx{0};
    for(const auto& at: step.getNames()){
        for(const auto& type: filter.types){
            if(at == type){
                tmp.push_back(idx);
                break;
            }
        }
        ++idx;
    }
    return tmp;
}

std::vector<size_t> evalIdx(const Step&, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    tmp = std::vector<size_t>{
            filter.indices.cbegin(), filter.indices.cend()};
    return tmp;
}

std::vector<size_t> evalPos(const Step& step, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    std::size_t idx{0};
    auto cmp = [&filter](const Vec& at){
        size_t dir = (filter.pos & filter.DIR_MASK) >> 3;
        if(filter.pos & filter.P_LT){
            return at[dir] < filter.posVal;
        }
        return at[dir] > filter.posVal;
    };
    for(const auto& at: step.asFmt(static_cast<AtomFmt>(filter.pos & filter.FMT_MASK))
                            .getCoords()){
        if(cmp(at)){
            tmp.push_back(idx);
        }
        ++idx;
    }
    return tmp;
}

std::vector<size_t> evalCoord(const Step& step, const SelectionFilter& filter){
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

std::vector<size_t> invertSel(const Step& step, const std::vector<size_t>& in){
    std::vector<size_t> out;
    out.reserve(step.getNat());
    for(size_t i=0; i<step.getNat(); ++i){
        if(std::find(in.begin(), in.end(), i) == in.end()){
            out.push_back(i);
        }
    }
    return out;
}

std::vector<size_t> evalFilter(const Step& step, SelectionFilter& filter);
std::vector<size_t> evalSubFilter(const Step& step,
                                  const SelectionFilter& filter,
                                  const std::vector<size_t>& parent){
    using Op = SelectionFilter::Op;
    std::vector<size_t> child = evalFilter(step, *filter.subfilter);
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

std::vector<size_t> evalFilter(const Step& step, SelectionFilter& filter)
{
    std::vector<size_t> tmp;
    switch(filter.mode){
    case SelectionFilter::Mode::Group:
        //TODO
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
    }
    if(filter.op & filter.NOT){
        tmp = invertSel(step, tmp);
    }
    if(filter.op & filter.PAIR){
        tmp = evalSubFilter(step, filter, tmp);
    }
    filter.op &= ~filter.UPDATE;
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
