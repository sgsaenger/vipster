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
    constexpr char fmtToC[] = {'b', 'a', 'c', 'd'};
    constexpr char dirToC[] = {'x', 'y', 'z'};
    os << dirToC[(filter.pos & Pos::DIR_MASK)>>2];
    os << fmtToC[filter.pos & Pos::FMT_MASK];
    os << ((filter.pos & Pos::P_CMP_MASK)? '<' : '>');
    os << filter.posVal;
}

static void writeIdx(std::ostream& os, const SelectionFilter& filter){
    os << "index ";
    if(filter.indices.size() == 1){
        os << filter.indices.cbegin()->first;
    }else{
        os << "[ ";
        for(const auto& p:filter.indices){
            os << p.first << " ";
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
        throw Error(std::string{"Invalid comparison "}+ctoken);
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
                filter.indices.emplace(static_cast<size_t>(i), std::vector{SizeVec{}});
            }
        }else{
            filter.indices.emplace(std::stoul(token), std::vector{SizeVec{}});
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
    filter = SelectionFilter{};
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
            throw Error("Invalid selection operator: "+token);
        }
    }
    // Parse subfilter (same level)
    if(is.good() && (is >> token) && !token.empty()){
        if(token.front() == ')'){
            is.putback(')');
            return is;
        }
        if(token.size() > 2){
            throw Error("Invalid coupling operator "+token);
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
            throw Error("Invalid coupling operator "+token);
        }
        filter.subfilter = std::make_unique<SelectionFilter>();
        is >> *filter.subfilter;
    }
    return is;
}
