#ifndef LIBVIPSTER_FILTER_TPP
#define LIBVIPSTER_FILTER_TPP

#include "filter.h"

using namespace Vipster;

template<typename T>
static SelectionIndices evalType(const T& step, const SelectionFilter& filter)
{
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
static SelectionIndices evalIdx(const T& step, const SelectionFilter& filter)
{
    SelectionIndices tmp;
    auto nat = step.getNat();
    for(const auto& p: filter.indices){
        // include all indices that refer to valid atoms
        if(p.first<nat){
            tmp.push_back(p);
        }
    }
    return tmp;
}

template<typename T>
static SelectionIndices evalPos(const T& step, const SelectionFilter& filter)
{
    // create comparison function
    std::function<bool(const Vec&)> cmp;
    size_t dir = (filter.pos & filter.DIR_MASK);
    switch(filter.pos & filter.P_CMP_MASK){
    case SelectionFilter::P_GT:
        cmp = [&](const Vec& at){return at[dir] > filter.posVal;};
        break;
    case SelectionFilter::P_LT:
        cmp = [&](const Vec& at){return at[dir] < filter.posVal;};
        break;
    case SelectionFilter::P_GEQ:
        cmp = [&](const Vec& at){return at[dir] >= filter.posVal;};
        break;
    case SelectionFilter::P_LEQ:
        cmp = [&](const Vec& at){return at[dir] <= filter.posVal;};
        break;
    default:
        //should be unreachable
        throw std::invalid_argument{"Invalid filter value"};
    }

    // evaluate criterion
    SelectionIndices tmp;
    std::size_t idx{0};
    for(const auto& at: step){
        if(cmp(at.coord)){
            tmp.emplace_back(idx, SizeVec{});
        }
        ++idx;
    }
    return tmp;
}

template<typename T>
static SelectionIndices evalCoord(const T& step, const SelectionFilter& filter)
{
    // determine coordination numbers
    // TODO: move functionality to step?
    std::vector<size_t> coord_numbers(step.getNat());
    for(const Bond& b: step.getBonds()){
        coord_numbers[b.at1] += 1;
        coord_numbers[b.at2] += 1;
    }

    // create comparison function
    std::function<bool(size_t)> cmp;
    auto cmp_op = filter.coord & filter.C_CMP_MASK;
    if(cmp_op == filter.C_GT){
        cmp = [&](size_t c){return c > filter.coordVal;};
    }else if(cmp_op == filter.C_EQ){
        cmp = [&](size_t c){return c == filter.coordVal;};
    }else{
        cmp = [&](size_t c){return c < filter.coordVal;};
    }

    // evaluate criterion
    SelectionIndices tmp;
    for(size_t i=0; i < step.getNat(); ++i){
        if(cmp(coord_numbers[i])){
            tmp.emplace_back(i, SizeVec{});
        }
    }

    return tmp;
}

template<typename T>
static SelectionIndices invertSel(const T& step, const SelectionIndices& in)
{
    SelectionIndices out{};
    const auto nat = step.getNat();
    out.reserve(nat);

    for(size_t i=0; i < nat; ++i){
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
                                  const SelectionIndices& parent)
{
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

#endif // LIBVIPSTER_FILTER_TPP
