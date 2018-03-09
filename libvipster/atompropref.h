#ifndef ATOMPROPREF_H
#define ATOMPROPREF_H

#include <iostream>
#include <utility>

namespace Vipster{

    // helper struct because std::bitset has no const_reference member type
    template<typename T>
    struct const_ref {
        using type = decltype(std::declval<const T>()[0]);
    };

    /*
     * Wrapper for Atom-properties
     *
     * dereferences pointers and exposes reference-semantics
     */
    class Atom;
    template<typename T, T* Atom::*p_prop, bool* Atom::*p_mod>
    class PropRef {
    public:
        PropRef(Atom& at):at{at} {}
        // can only be explicitely constructed to reference an Atom
        PropRef(const PropRef&) = delete;
        // like a real reference, assigning changes the origin
        PropRef& operator=(const PropRef& rhs){
            *(at.*p_prop) = static_cast<T>(rhs);
            *(at.*p_mod) = true;
            return *this;
        }
        // convert to reference
        operator T&() {*(at.*p_mod) = true; return *(at.*p_prop);}
        operator const T&() const {return *(at.*p_prop);}
        PropRef& operator=(const T& prop)
        {
            *(at.*p_prop) = prop;
            *(at.*p_mod) = true;
            return *this;
        }
        // enable c_str-conversion for name
        template <typename R=T>
        const typename R::value_type* c_str() const {return (at.*p_prop)->c_str();}
        // enable array-access for Vec
        template <typename R=T>
        typename R::reference operator [](size_t i) {
            *(at.*p_mod) = true;
            return (*(at.*p_prop))[i];
        }
        template <typename R=T>
        typename const_ref<R>::type operator [](size_t i) const {
            return (*(at.*p_prop))[i];
        }
        bool operator==(const PropRef &rhs) const
        {
            return *(at.*p_prop) == static_cast<T>(rhs);
        }
    private:
        Atom& at;
    };
    template<typename T, T* Atom::*p_prop, bool* Atom::*p_mod>
    std::ostream& operator<<(std::ostream& s, const PropRef<T,p_prop,p_mod>& pr){
        s << static_cast<const T&>(pr);
        return s;
    }
    template<typename T, T* Atom::*p_prop, bool* Atom::*p_mod>
    std::istream& operator>>(std::istream& s, PropRef<T,p_prop,p_mod>& pr){
        s >> static_cast<T&>(pr);
        return s;
    }
    template<typename T, T* Atom::*p_prop, bool* Atom::*p_mod>
    bool operator==(const T& lhs, const PropRef<T,p_prop,p_mod>& rhs){
        return lhs == static_cast<const T&>(rhs);
    }
    template<typename T, T* Atom::*p_prop, bool* Atom::*p_mod>
    bool operator==(const PropRef<T,p_prop,p_mod>& lhs, const T& rhs){
        return static_cast<const T&>(lhs) == rhs;
    }

}

#endif
