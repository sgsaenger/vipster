#ifndef STEPFMT_H
#define STEPFMT_H

#include "atom.h"

#include <functional>
#include <memory>

namespace Vipster::detail{

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

    AtomContext ctxt;
    // store original atom-storage
    std::shared_ptr<T> atoms;
    // complete Step-interface
    size_t getNat() const {return atoms->getNat();}
    // store format functions
    fmtfun ffun, invfun;

    Formatter(const std::shared_ptr<T> &atoms, AtomFmt fmt)
        : atoms{atoms},
          ctxt{fmt, atoms->ctxt.pte, atoms->ctxt.cell},
          ffun{makeConverter(atoms->ctxt, ctxt)},
          invfun{makeConverter(ctxt, atoms->ctxt)}
    {}
    Formatter(const Formatter&) = default;
    Formatter(Formatter&&) = default;
    Formatter& operator=(const Formatter&) = default;
    Formatter& operator=(Formatter&&) = default;

    template<bool isConst>
    class AtomView : private T::template AtomView<isConst>{
        /* Wrapper class that wraps existing Atom accessor
         * and transforms coordinate accesses as needed
         */
        template<bool> friend class AtomView;
        using base = typename T::template AtomView<isConst>;
        AtomView()
            : base{}, source{nullptr}
        {}
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
            _Vec& operator=(const _Vec& rhs){
                return operator=(static_cast<const Vec&>(rhs));
            }
            _Vec& operator=(const Vec& rhs){
                a.base::coord = a.source->invfun(rhs);
                return *this;
            }
            // implement modify-assignment here because we can't convert to reference
            _Vec& operator+=(const Vec& rhs){
                a.base::coord += a.source->invfun(rhs);
                return *this;
            }
            _Vec& operator-=(const Vec& rhs){
                a.base::coord -= a.source->invfun(rhs);
                return *this;
            }
            _Vec& operator+=(const double &f){
                a.base::coord += a.source->invfun({f,f,f});
                return *this;
            }
            _Vec& operator-=(const double &f){
                a.base::coord += a.source->invfun({f,f,f});
                return *this;
            }
            _Vec& operator*=(const double &f){
                a.base::coord *= f;
                return *this;
            }
            _Vec& operator/=(const double &f){
                a.base::coord /= f;
                return *this;
            }
            // convert to Vec
            operator Vec() const {return a.source->ffun(a.base::coord);}
            Vec asFmt(const AtomContext &ctxt) const
            {
                return makeConverter(a.source->ctxt, ctxt)(static_cast<const Vec&>(*this));
            }
            /* array access
             * const-access as usual,
             * mutable access to another helper that wraps assignment between back-and-forth conversion
             */
            const Vec::value_type& operator[](std::size_t i) const{
                return a.source->ffun(a.base::coord)[i];
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
            bool operator==(const Vec &rhs) const {return rhs == operator Vec();}
        };
    protected:
        Formatter *source;
        AtomView& operator+=(ptrdiff_t i){
            base::operator+=(i);
            return *this;
        }
        void pointTo(const AtomView &rhs){
            base::pointTo(rhs);
            source = rhs.source;
        }
    public:
        // "Data", encapsulated in wrapper objects
        std::conditional_t<isConst, const _Vec, _Vec>   coord{*this};
        using base::name;
        using base::type;
        using base::properties;
        using base::operator==;

        AtomView(Formatter &f, size_t i)
            : base{*f.atoms, i}, source{&f}
        {}
        // copying is templated to allow conversion to const
        // copy constructor creates new object pointing to same data
        AtomView(const AtomView &rhs)
            : base{rhs}, source{rhs.source}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomView(const AtomView &rhs)
            : base{rhs}, source{rhs.source}
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
        virtual ~AtomView() = default;
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
            : value_type{}, idx{}
        {}
        AtomIterator(Formatter &f, size_t i)
            : value_type{f, i}, idx{i}
        {}
        // copying is templated to allow conversion to const
        // default copy constructor
        AtomIterator(const AtomIterator &it)
            : value_type{it}, idx{it.idx}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomIterator(const AtomIterator<B> &it)
            : value_type{it}, idx{it.idx}
        {}
        // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomIterator&   operator=(const AtomIterator &it){
            value_type::pointTo(it);
            idx = it.idx;
            return *this;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
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
    };
};

}

#endif // STEPFMT_H
