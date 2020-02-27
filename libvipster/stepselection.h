#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include "atom.h"
#include "filter.h"
#include <vector>
#include <memory>

namespace Vipster {
/*
 * Selection container
 *
 * contains indices of selected atoms in Step-like
 */
template<typename T>
struct Selection{
    template<bool isConst>
    class AtomView;
    using atom = AtomView<false>;
    using const_atom = AtomView<true>;
    template<bool isConst>
    class AtomIterator;
    using iterator = AtomIterator<false>;
    using const_iterator = AtomIterator<true>;

    Selection(std::shared_ptr<T> atoms, SelectionIndices indices)
        : atoms{atoms}, fmt{atoms->fmt}, indices{indices}
    {}
    Selection(const Selection&) = default;
    Selection(Selection &&) = default;
    Selection& operator=(const Selection&) = default;
    Selection& operator=(Selection&&) = default;

    // store original atom-storage
    std::shared_ptr<T> atoms;
    // complete Step-interface
    AtomFmt fmt;
    size_t getNat();
    // store selection
    SelectionIndices indices;

    template<bool isConst>
    class AtomView: protected T::template AtomView<isConst>
    {
        using base = typename T::template AtomView<isConst>;
    public:
        AtomView(Selection &s, PeriodicTable &pte, size_t i)
            : base{*s.atoms, pte, i}
        {}
        using base::coord;
        using base::name;
        using base::type;
        using base::properties;
    };

    template<bool isConst>
    class AtomIterator: private AtomView<isConst>
    {
    public:
        using difference_type = ptrdiff_t;
        using value_type = AtomView<isConst>;
        using reference = value_type&;
        using pointer = value_type*;
        using iterator_category = std::random_access_iterator_tag;
        // TODO: default constructibility
        AtomIterator(Selection &s, PeriodicTable &pte, size_t i)
            : value_type{s, pte, s.indices[i].first},
              idx{i}, sel{&s.indices}
        {}
        // copying is templated to allow conversion to const
        // default copy constructor
        AtomIterator(const AtomIterator &it)
            : value_type{it}, idx{it.idx}, sel{it.sel}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t||!B>::type>
        AtomIterator(const AtomIterator<B> &it)
            : value_type{it}, idx{it.idx}, sel{it.sel}
        {}
        // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomIterator&   operator=(const AtomIterator &it){
            idx = it.idx;
            sel = it.sel;
            this->v = it.v;
            this->elem = it.elem;
            this->pte = it.pte;
            this->prop = it.prop;
            return *this;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t||!B>::type>
        AtomIterator&   operator=(const AtomIterator<B> &it){
            idx = it.idx;
            sel = it.sel;
            this->v = it.v;
            this->elem = it.elem;
            this->pte = it.pte;
            this->prop = it.prop;
            return *this;
        }
        // access
        reference   operator*() const {
            // remove constness of iterator, as it is independent of constness of Atom
            return static_cast<reference>(*const_cast<AtomIterator*>(this));
        }
        pointer     operator->() const {
            return &(operator*());
        }
        reference   operator[](difference_type i){
            return *operator+(i);
        }
        size_t      getIdx() const noexcept{
            return idx;
        }
        const SelectionPair& getSel() const{
            return (*sel)[idx];
        }
        // comparison
        difference_type operator-(const AtomIterator &rhs) const{
            return getIdx() - rhs.getIdx();
        }
        bool            operator==(const AtomIterator &rhs) const{
            return (this->elem == rhs.elem);
        }
        bool            operator!=(const AtomIterator &rhs) const{
            return !operator==(rhs);
        }
        friend bool     operator< (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx < rhs.idx;
        }
        friend bool     operator> (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx > rhs.idx;
        }
        friend bool     operator<= (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx <= rhs.idx;
        }
        friend bool     operator>= (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx >= rhs.idx;
        }
        // in-/decrement
        AtomIterator&   operator++(){
            operator+=(1);
            return *this;
        }
        AtomIterator    operator++(int){
            auto copy = *this;
            operator+=(1);
            return copy;
        }
        AtomIterator&   operator+=(difference_type i){
            idx += i;
            auto diff = (*sel)[i].first - (*sel)[idx].first;
            this->v += diff;
            this->elem += diff;
            this->prop += diff;
            return *this;
        }
        AtomIterator    operator+(difference_type i){
            auto copy = *this;
            return copy+=i;
        }
        AtomIterator&   operator--(){
            operator+=(-1);
            return *this;
        }
        AtomIterator    operator--(int){
            auto copy = *this;
            operator+=(-1);
            return copy;
        }
        AtomIterator&   operator-=(difference_type i){
            operator+=(-i);
            return *this;
        }
        AtomIterator    operator-(difference_type i){
            auto copy = *this;
            return copy-=i;
        }
    private:
        size_t idx;
        std::vector<SelectionPair> *sel;
    };
};

}

#endif // LIBVIPSTER_STEPSEL_H
