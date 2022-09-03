#ifndef ATOM_TPP
#define ATOM_TPP

#include <cstddef>
#include <iterator>

#include "atom.h"

using namespace Vipster::detail;

// Construction
template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>::AtomIterator()
    : value_type{}, idx{}
{}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>::AtomIterator(typename value_type::Source &s, size_t i)
    : value_type{s, i}, idx{i}
{}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>::AtomIterator(const AtomIterator &it)
    : value_type{it}, idx{it.idx}
{}

template<template<bool> typename AtomView, bool isConst>
template<bool B, bool t, typename>
AtomIterator<AtomView, isConst>::AtomIterator(const AtomIterator<AtomView, B> &it)
    : value_type{it}, idx{it.idx}
{}

// Assignment
template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>& AtomIterator<AtomView, isConst>::operator=(const AtomIterator &it){
    value_type::pointTo(it);
    idx = it.idx;
    return *this;
}

template<template<bool> typename AtomView, bool isConst>
template<bool B, bool t, typename>
AtomIterator<AtomView, isConst>& AtomIterator<AtomView, isConst>::operator=(const AtomIterator<AtomView, B> &it){
    value_type::pointTo(it);
    idx = it.idx;
    return *this;
}

// access
template<template<bool> typename AtomView, bool isConst>
typename AtomIterator<AtomView, isConst>::reference   AtomIterator<AtomView, isConst>::operator*() const {
    // remove constness of iterator, as it is independent of constness of Atom
    return static_cast<reference>(*const_cast<AtomIterator*>(this));
}

template<template<bool> typename AtomView, bool isConst>
typename AtomIterator<AtomView, isConst>::pointer     AtomIterator<AtomView, isConst>::operator->() const {
    return &(operator*());
}

template<template<bool> typename AtomView, bool isConst>
typename AtomIterator<AtomView, isConst>::reference   AtomIterator<AtomView, isConst>::operator[](difference_type i){
    return *operator+(i);
}

template<template<bool> typename AtomView, bool isConst>
size_t      AtomIterator<AtomView, isConst>::getIdx() const noexcept{
    return idx;
}

// comparison
template<template<bool> typename AtomView, bool isConst>
typename AtomIterator<AtomView, isConst>::difference_type AtomIterator<AtomView, isConst>::operator-(const AtomIterator &rhs) const{
    return getIdx() - rhs.getIdx();
}

template<template<bool> typename AtomView, bool isConst>
bool            AtomIterator<AtomView, isConst>::operator==(const AtomIterator &rhs) const{
    return (this->source == rhs.source) && (this->idx == rhs.idx);
}

template<template<bool> typename AtomView, bool isConst>
bool            AtomIterator<AtomView, isConst>::operator!=(const AtomIterator &rhs) const{
    return !operator==(rhs);
}

template<template<bool> typename AtomView, bool isConst>
bool            AtomIterator<AtomView, isConst>::operator< (const AtomIterator &rhs) const{
    return idx < rhs.idx;
}

template<template<bool> typename AtomView, bool isConst>
bool            AtomIterator<AtomView, isConst>::operator> (const AtomIterator &rhs) const{
    return idx > rhs.idx;
}

template<template<bool> typename AtomView, bool isConst>
bool            AtomIterator<AtomView, isConst>::operator<=(const AtomIterator &rhs) const{
    return idx <= rhs.idx;
}

template<template<bool> typename AtomView, bool isConst>
bool            AtomIterator<AtomView, isConst>::operator>=(const AtomIterator &rhs) const{
    return idx >= rhs.idx;
}

// in-/decrement
template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>&   AtomIterator<AtomView, isConst>::operator++(){
    operator+=(1);
    return *this;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>    AtomIterator<AtomView, isConst>::operator++(int){
    auto copy = *this;
    operator+=(1);
    return copy;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>&   AtomIterator<AtomView, isConst>::operator+=(difference_type i){
    idx += i;
    value_type::operator+=(i);
    return *this;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>    AtomIterator<AtomView, isConst>::operator+(difference_type i){
    auto copy = *this;
    return copy+=i;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>&   AtomIterator<AtomView, isConst>::operator--(){
    operator+=(-1);
    return *this;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>    AtomIterator<AtomView, isConst>::operator--(int){
    auto copy = *this;
    operator+=(-1);
    return copy;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst>&   AtomIterator<AtomView, isConst>::operator-=(difference_type i){
    operator+=(-i);
    return *this;
}

template<template<bool> typename AtomView, bool isConst>
AtomIterator<AtomView, isConst> AtomIterator<AtomView, isConst>:: operator-(difference_type i){
    auto copy = *this;
    return copy-=i;
}

#endif // ATOM_TPP
