#ifndef STEPFMT_H
#define STEPFMT_H

#include "atom.h"

#include <memory>

namespace Vipster{

template <typename T>
struct Formatter{
    template<bool isConst>
    class AtomView;
    using atom = AtomView<false>;
    using const_atom = AtomView<true>;
    template<bool isConst>
    class AtomIterator;
    using iterator = AtomIterator<false>;
    using const_iterator = AtomIterator<true>;
    using fmtfun = std::function<Vec(const Vec&)>;

    Formatter(const std::shared_ptr<T> &atoms, AtomFmt fmt, fmtfun ffun, fmtfun inv)
        : atoms{atoms}, fmt{fmt}, ffun{ffun}, invfun{inv}
    {}
    Formatter(const Formatter&) = default;
    Formatter(Formatter&&) = default;
    Formatter& operator=(const Formatter&) = default;
    Formatter& operator=(Formatter&&) = default;

    // store original atom-storage
    std::shared_ptr<T> atoms;
    // complete Step-interface
    AtomFmt fmt;
    size_t getNat() const {return atoms->getNat();}
    // store format functions
    fmtfun ffun, invfun;

    template<bool isConst>
    class AtomView : protected T::template AtomView<isConst>{
        /* Wrapper class that wraps existing Atom accessor
         * and transforms coordinate accesses as needed
         */
        using base = typename T::template AtomView<isConst>;
        class _Vec{
        private:
            AtomView &a;
        public:
            _Vec(AtomView &a)
                : a{a}
            {}
            // only explicit initialization
            _Vec(const _Vec&) = delete;
            // assigning changes the origin
            _Vec& operator=(const Vec& rhs){
                *a.v = a.inv(rhs);
                return *this;
            }
            // implement modify-assignment here because we can't convert to reference
            _Vec& operator+=(const Vec& rhs){
                *a.v += a.inv(rhs);
                return *this;
            }
            _Vec& operator-=(const Vec& rhs){
                *a.v -= a.inv(rhs);
                return *this;
            }
            _Vec& operator+=(const double &f){
                *a.v += a.inv({f,f,f});
                return *this;
            }
            _Vec& operator-=(const double &f){
                *a.v += a.inv({f,f,f});
                return *this;
            }
            _Vec& operator*=(const double &f){
                *a.v *= f;
                return *this;
            }
            _Vec& operator/=(const double &f){
                *a.v /= f;
                return *this;
            }
            // convert to Vec
            operator Vec() const {return a.fmt(*a.v);}
            /* array access
             * const-access as usual,
             * mutable access to another wrapper that stores the formatted vec and reassigns to the origin
             */
            const Vec::value_type& operator[](std::size_t i) const{
                return a.fmt(*a.v)[i];
            }
            class _Value_type{
            public:
                _Value_type(_Vec& _v, size_t i): v{_v.a.fmt(_v)}, _v{_v}, i{i} {}
                operator const Vec::value_type&() const {return v[i];}
                _Value_type& operator=(Vec::value_type val){
                    v[i] = val;
                    _v = _v.a.inv(v);
                    v = _v.a.fmt(_v); // reload the local vector to minimize race-conditions if multiple _Value_types exist
                    return *this;
                }
                _Value_type& operator-=(Vec::value_type val){
                    v[i] -= val;
                    _v = _v.a.inv(v);
                    v = _v.a.fmt(_v); // reload the local vector to minimize race-conditions if multiple _Value_types exist
                    return *this;
                }
                _Value_type& operator+=(Vec::value_type val){
                    v[i] += val;
                    _v = _v.a.inv(v);
                    v = _v.a.fmt(_v); // reload the local vector to minimize race-conditions if multiple _Value_types exist
                    return *this;
                }
            private:
                Vec v;
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
        fmtfun &fmt, &inv;
    public:
        // "Data", encapsulated in wrapper objects
        std::conditional_t<isConst, const _Vec, _Vec>   coord{*this};
        using base::name;
        using base::type;
        using base::properties;

        AtomView(Formatter &f, PeriodicTable &pte, size_t i)
            : base{*f.atoms, pte, i}, fmt{f.ffun}, inv{f.invfun}
        {}
        AtomView(const AtomView &rhs)
            : base{rhs}, fmt{rhs.fmt}, inv{rhs.inv}
        {}
        virtual ~AtomView() = default;
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
        AtomIterator(Formatter &f, PeriodicTable &pte, size_t i)
            : value_type{f, pte, i}, idx{i}
        {}
        // copying is templated to allow conversion to const
        // default copy constructor
        AtomIterator(const AtomIterator &it)
            : value_type{it}, idx{it.idx}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t||!B>::type>
        AtomIterator(const AtomIterator<B> &it)
            : value_type{it}, idx{it.idx}
        {}
        // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomIterator&   operator=(const AtomIterator &it){
            this->v = it.v;
            this->elem = it.elem;
            this->pte = it.pte;
            this->prop = it.prop;
            idx = it.idx;
            return *this;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t||!B>::type>
        AtomIterator&   operator=(const AtomIterator<B> &it){
            this->v = it.v;
            this->elem = it.elem;
            this->pte = it.pte;
            this->prop = it.prop;
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
            this->v += i;
            this->elem += i;
            this->prop += i;
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
    };
};

}

#endif // STEPFMT_H
