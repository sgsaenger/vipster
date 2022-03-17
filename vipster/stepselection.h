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
    using AtomIterator = detail::AtomIterator<AtomView, isConst>;
    using iterator = AtomIterator<false>;
    using const_iterator = AtomIterator<true>;

    AtomContext ctxt, crys_ctxt;
    std::function<Vec(const Vec&)> ffun;
    // store original atom-storage
    std::shared_ptr<T> atoms;
    // complete Step-interface
    size_t getNat() const noexcept {return indices.size();}
    // store selection
    SelectionIndices indices;

    Selection(std::shared_ptr<T> atoms, SelectionIndices indices)
        : ctxt{atoms->ctxt},
          atoms{atoms},
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
};

template<typename T>
template<bool isConst>
class Selection<T>::AtomView: private T::template AtomView<isConst>
{
    SelectionPair *sel;
    template<bool> friend class AtomView;
    using base = typename T::template AtomView<isConst>;
    AtomView()
        : base{}, source{}, sel{}
    {}
    class _Index{
    private:
        AtomView &a;
    public:
        _Index(AtomView &a)
            : a{a}
        {}
        // only explicit initialization, no assignment
        _Index(const _Index&) = delete;
        // expose index in source-Step
        operator int() { return a.sel->first; }
    };
    class _Offset{
    private:
        AtomView &a;
    public:
        _Offset(AtomView &a)
            : a{a}
        {}
        // only explicit initialization, no assignment
        _Offset(const _Offset&) = delete;
        // assigning changes the origin
        _Offset& operator=(const _Offset &rhs){
            a.sel->second = rhs;
            return *this;
        }
        _Offset& operator=(const SizeVec &rhs){
            a.sel->second = rhs;
            return *this;
        }
        // Convert to SizeVec -> cell-based offset
        operator const SizeVec&() const {
            return a.sel->second;
        }
        // Convert to Vec -> real offset in source's format
        operator Vec() const {
            const auto& off = a.sel->second;
            return a.source->ffun(Vec{static_cast<Vec::value_type>(off[0]),
                                      static_cast<Vec::value_type>(off[1]),
                                      static_cast<Vec::value_type>(off[2])});
        }
        // array access
        const SizeVec::value_type& operator[](std::size_t i) const {
            return a.sel->second[i];
        }
        SizeVec::value_type& operator[](std::size_t i){
            return a.sel->second[i];
        }
    };
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
            a.base::coord = rhs - a.off;
            return *this;
        }
        // implement modify-assignment here because we can't convert to reference
        _Vec& operator+=(const Vec& rhs){
            a.base::coord += rhs;
            return *this;
        }
        _Vec& operator-=(const Vec& rhs){
            a.base::coord -= rhs;
            return *this;
        }
        _Vec& operator+=(const double &f){
            a.base::coord += Vec{f,f,f};
            return *this;
        }
        _Vec& operator-=(const double &f){
            a.base::coord -= Vec{f,f,f};
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
        operator Vec() const {return a.base::coord + a.off;}
        Vec asFmt(const AtomContext &ctxt) const
        {
            return makeConverter(a.source->ctxt, ctxt)(*this);
        }
        /* array access
         * const-access as usual,
         * mutable access to another helper that wraps assignment between back-and-forth conversion
         */
        Vec::value_type operator[](std::size_t i) const{
            return (a.base::coord + a.off)[i];
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
    using Source = Selection;
    Source *source;
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
    _Offset off{*this};
    _Index idx{*this};
    std::conditional_t<isConst, const _Vec, _Vec>   coord{*this};
    using base::name;
    using base::type;
    using base::properties;
    using base::operator==;

    AtomView(Selection &s, size_t i)
        : base{*s.atoms, i>=s.indices.size() ? 0 : s.indices[i].first},
          sel{&*(s.indices.begin()+i)},
          source{&s}
    {}
    // copying is templated to allow conversion to const
    // copy constructor creates new object pointing to same data
    AtomView(const AtomView &rhs)
        : base{rhs}, sel{rhs.sel}, source{rhs.source}
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

}

#endif // LIBVIPSTER_STEPSEL_H
