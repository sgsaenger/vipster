#include "stepsel.h"
#include "step.h"

using namespace Vipster;

static void writeType(std::ostream& os, const SelectionFilter& filter){
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

static void writeCoord(std::ostream& os, const SelectionFilter& filter){
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

static void writePos(std::ostream& os, const SelectionFilter& filter){
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

static void writeIdx(std::ostream& os, const SelectionFilter& filter){
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

std::ostream& Vipster::operator<<(std::ostream& os, const SelectionFilter& filter)
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
    default:
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

static void parseType(std::istream& is, SelectionFilter& filter){
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

static void parseCoord(std::istream& is, SelectionFilter& filter){
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

static void parsePos(std::istream& is, SelectionFilter& filter){
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

static void parseIdx(std::istream& is, SelectionFilter& filter){
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

std::istream& Vipster::operator>>(std::istream& is, SelectionFilter& filter){
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
static std::vector<size_t> evalType(const Step& step, const SelectionFilter& filter){
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

static std::vector<size_t> evalIdx(const Step&, const SelectionFilter& filter){
    std::vector<size_t> tmp;
    tmp = std::vector<size_t>{
            filter.indices.cbegin(), filter.indices.cend()};
    return tmp;
}

static std::vector<size_t> evalPos(const Step& step, const SelectionFilter& filter){
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

static std::vector<size_t> evalCoord(const Step& step, const SelectionFilter& filter){
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


static std::vector<size_t> invertSel(const Step& step, const std::vector<size_t>& in){
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

static std::vector<size_t> evalSubFilter(const Step& step,
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

std::vector<size_t> Vipster::evalFilter(const Step& step, SelectionFilter& filter)
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
