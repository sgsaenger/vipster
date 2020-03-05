#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include "atom.h"
#include "filter.h"

#include <memory>

namespace Vipster::detail{

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
    using fmtfun = std::function<Vec(const Vec&)>;

    AtomContext ctxt, crys_ctxt;
    fmtfun ffun;
    // store original atom-storage
    std::shared_ptr<T> atoms;
    // complete Step-interface
    size_t getNat() const noexcept {return indices.size();}
    // store selection
    SelectionIndices indices;

    Selection(std::shared_ptr<T> atoms, SelectionIndices indices)
        : atoms{atoms},
          ctxt{atoms->ctxt},
          indices{indices}
    {
        crys_ctxt = ctxt;
        crys_ctxt.fmt = AtomFmt::Crystal;
        ffun = makeConverter(crys_ctxt, ctxt);
    }
    Selection(const Selection&) = default;
    Selection(Selection &&) = default;
    Selection& operator=(const Selection&) = default;
    Selection& operator=(Selection&&) = default;

    template<bool isConst>
    class AtomView: private T::template AtomView<isConst>
    {
        SelectionIndices::value_type *sel;
        template<bool> friend class AtomView;
        using base = typename T::template AtomView<isConst>;
        AtomView()
            : base{}, source{}, sel{}
        {}
        class _Vec{
        private:
            AtomView &a;
            Vec makeOff() const {
                const auto& off = a.sel->second;
                return a.source->ffun(Vec{static_cast<Vec::value_type>(off[0]),
                                          static_cast<Vec::value_type>(off[1]),
                                          static_cast<Vec::value_type>(off[2])});
            }
        public:
            _Vec(AtomView &a)
                : a{a}
            {}
            // only explicit initialization
            _Vec(const _Vec&) = delete;
            // assigning changes the origin
            _Vec& operator=(const _Vec& rhs){
                return operator=(static_cast<const Vec&>(rhs));
            }
            _Vec& operator=(const Vec& rhs){
                a.base::coord = rhs - makeOff();
                return *this;
            }
            // implement modify-assignment here because we can't convert to reference
            _Vec& operator+=(const Vec& rhs){
                a.base::coord += rhs - makeOff();
                return *this;
            }
            _Vec& operator-=(const Vec& rhs){
                a.base::coord -= rhs - makeOff();
                return *this;
            }
            _Vec& operator+=(const double &f){
                a.base::coord += Vec{f,f,f} - makeOff();
                return *this;
            }
            _Vec& operator-=(const double &f){
                a.base::coord -= Vec{f,f,f} - makeOff();
                return *this;
            }
            _Vec& operator*=(const double &f){
                a.base::coord *= f;
                return *this;
            }
            _Vec& operator/=(const double &f){
                a.base::coord *= f;
                return *this;
            }
            // convert to Vec
            operator Vec() const {return a.base::coord + makeOff();}
            Vec asFmt(const AtomContext &ctxt) const
            {
                return makeConverter(*a.ctxt, ctxt)(static_cast<const Vec&>(*this));
            }
            /* array access
             * const-access as usual,
             * mutable access to another helper that wraps assignment between back-and-forth conversion
             */
            const Vec::value_type& operator[](std::size_t i) const{
                return (a.base::coord + makeOff())[i];
            }
            class _Value_type{
            public:
                _Value_type(_Vec& _v, size_t i): _v{_v}, i{i} {}
                operator const Vec::value_type&() const {return _v[i];}
                _Value_type& operator=(Vec::value_type val){
                    Vec v = _v;
                    v[i] = val;
                    _v = v;
                    return *this;
                }
                _Value_type& operator-=(Vec::value_type val){
                    Vec v = _v;
                    v[i] -= val;
                    _v = v;
                    return *this;
                }
                _Value_type& operator+=(Vec::value_type val){
                    Vec v = _v;
                    v[i] += val;
                    _v = v;
                    return *this;
                }
            private:
                _Vec &_v;
                size_t i;
            };
            _Value_type operator[] (std::size_t i){
                return {*this, i};
            }
            // comparison
            bool operator==(const Vec &rhs) const {return rhs == *this;}
        };
    protected:
        Selection *source;
        AtomView& operator+=(ptrdiff_t i){
            auto diff = (sel+i)->first - sel->first;
            base::operator+=(diff);
            this->sel += i;
            return *this;
        }
        void pointTo(const AtomView &rhs){
            base::pointTo(rhs);
            sel = rhs.sel;
            source = rhs.source;
        }
    public:
        // "Data", encapsulated in wrapper objects
        std::conditional_t<isConst, const _Vec, _Vec>   coord{*this};
        using base::name;
        using base::type;
        using base::properties;
        using base::operator==;

        AtomView(Selection &s, size_t i)
            : base{*s.atoms, s.indices[i].first}, source{&s}, sel{&s.indices[i]}
        {}
        // copying is templated to allow conversion to const
        // copy constructor creates new object pointing to same data
        AtomView(const AtomView &rhs)
            : base{rhs}, source{rhs.source}, sel{rhs.sel}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomView(const AtomView &rhs)
            : base{rhs}, source{rhs.source}, sel{rhs.sel}
        {}
        // copy assignment changes data
        AtomView& operator=(const AtomView &rhs){
            coord = rhs.coord.asFmt(source->ctxt);
            name = rhs.name;
            properties = rhs.properties;
            return *this;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomView& operator=(const AtomView &rhs){
            coord = rhs.coord.asFmt(source->ctxt);
            name = rhs.name;
            properties = rhs.properties;
            return *this;
        }
        template<bool B, template<bool> typename AV>
        AtomView& operator=(const AV<B> &rhs){
            coord = rhs.coord.asFmt(source->ctxt);
            name = rhs.name;
            properties = rhs.properties;
            return *this;
        }
    };

    template<bool isConst>
    class AtomIterator: private AtomView<isConst>
    {
        template<bool> friend class AtomIterator;
    public:
        using difference_type = ptrdiff_t;
        using value_type = AtomView<isConst>;
        using reference = value_type&;
        using pointer = value_type*;
        using iterator_category = std::random_access_iterator_tag;

        AtomIterator()
            : value_type{}, idx{}, source{nullptr}
        {}
        AtomIterator(Selection &s, size_t i)
            : value_type{s, (i>=s.indices.size()) ? s.atoms->getNat() : s.indices[i].first},
              idx{i}, source{&s.indices}
        {}
        // copying is templated to allow conversion to const
        // default copy constructor
        AtomIterator(const AtomIterator &it)
            : value_type{it}, idx{it.idx}, source{it.source}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t||!B>::type>
        AtomIterator(const AtomIterator<B> &it)
            : value_type{it}, idx{it.idx}, source{it.source}
        {}
        // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomIterator&   operator=(const AtomIterator &it){
            value_type::pointTo(it);
            idx = it.idx;
            return *this;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t||!B>::type>
        AtomIterator&   operator=(const AtomIterator<B> &it){
            value_type::pointTo(it);
            idx = it.idx;
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
            return (*source)[idx];
        }
        // comparison
        difference_type operator-(const AtomIterator &rhs) const{
            return getIdx() - rhs.getIdx();
        }
        bool            operator==(const AtomIterator &rhs) const{
            return (this->source == rhs.source) && (this->idx == rhs.idx);
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
        friend bool     operator<=(const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx <= rhs.idx;
        }
        friend bool     operator>=(const AtomIterator &lhs, const AtomIterator &rhs){
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
            value_type::operator+=(i);
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
        std::vector<SelectionPair> *source;
    };
};

}

#endif // LIBVIPSTER_STEPSEL_H
