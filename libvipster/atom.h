#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <tuple>
#include <bitset>
#include <utility>
#include <iostream>

#include "vec.h"
#include "pse.h"

//TODO: track changes in properties?

namespace Vipster{
    namespace Enums{
    enum AtomChange{coord=0x1, name=0x2, prop=0x4, del=0x8};
    enum AtomFlag: uint8_t {FixX, FixY, FixZ, Hidden};
    }
    using Enums::AtomChange;
    using Enums::AtomFlag;
    constexpr size_t nAtFlag = 4;
    using AtomFlags = std::bitset<nAtFlag>;

    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };
    constexpr size_t nAtFmt = 4;

    struct AtomProperties{
        float charge;
        Vec forces;
        AtomFlags flags;
    };
    inline bool operator==(const AtomProperties &p1, const AtomProperties &p2){
        return std::tie(p1.charge, p1.flags, p1.forces)
               ==
               std::tie(p2.charge, p2.flags, p2.forces);
    }


    /*
     * Basic atom interface
     *
     * Provides access-wrappers to properties
     */
    template<typename T>
    class AtomViewBase{
        template<typename U> friend class AtomViewBase;
    protected:
        // used as base for iterators, no standalone creation intended atm.
        AtomViewBase(Vec *co, bool *c_m, std::string *n, bool *n_m,
             AtomProperties *p, PseEntry** pse, bool *p_m)
            : coord_ptr{co}, coord_changed{c_m}, name_ptr{n}, name_changed{n_m},
              prop_ptr{p}, pse_ptr{pse}, prop_changed{p_m}
        {}
        // allow Atom to constAtom conversion
        template <typename U, typename R=T, typename = typename std::enable_if<std::is_const<R>::value>::type>
        AtomViewBase(const AtomViewBase<U>& at)
            : coord_ptr{at.coord_ptr}, coord_changed{at.coord_changed},
              name_ptr{at.name_ptr}, name_changed{at.name_changed},
              prop_ptr{at.prop_ptr}, pse_ptr{at.pse_ptr}, prop_changed{at.prop_changed}
        {}
        Vec                     *coord_ptr;
        bool                    *coord_changed;
        std::string             *name_ptr;
        bool                    *name_changed;
        AtomProperties          *prop_ptr;
        PseEntry*               *pse_ptr;
        bool                    *prop_changed;
    private:
        template<typename U, U* AtomViewBase::*p_prop, bool* AtomViewBase::*p_mod>
        class PropRef;
        template<typename U, U* AtomViewBase::*p_prop, bool* AtomViewBase::*p_mod>
        using ref = typename std::conditional<std::is_const<T>::value,
                        const class PropRef<U,p_prop,p_mod>,
                        class PropRef<U,p_prop,p_mod>>::type;
    public:
        ref<std::string, &AtomViewBase::name_ptr, &AtomViewBase::name_changed>
            name{*this};
        ref<Vec, &AtomViewBase::coord_ptr, &AtomViewBase::coord_changed>
            coord{*this};
        ref<AtomProperties, &AtomViewBase::prop_ptr, &AtomViewBase::prop_changed>
            properties{*this};
        const PropRef<PseEntry*, &AtomViewBase::pse_ptr, nullptr>
            pse{*this};
        virtual ~AtomViewBase() = default;
        AtomViewBase() = delete;
        inline AtomViewBase(const AtomViewBase& at)
            : coord_ptr{at.coord_ptr},
              coord_changed{at.coord_changed},
              name_ptr{at.name_ptr},
              name_changed{at.name_changed},
              prop_ptr{at.prop_ptr},
              pse_ptr{at.pse_ptr},
              prop_changed{at.prop_changed} {}
        inline AtomViewBase& operator=(const AtomViewBase& at)
        {
            name = at.name;
            coord = at.coord;
            properties = at.properties;
            pse_ptr = at.pse_ptr;
            return *this;
        }
        template <typename U>
        bool operator==(const AtomViewBase<U>& rhs) const
        {
            return std::tie(name, coord, properties)
                   ==
                   std::tie(rhs.name, rhs.coord, rhs.properties);
        }
        template <typename U>
        bool operator!=(const AtomViewBase<U>& rhs) const
        {
            return !(*this==rhs);
        }
    private:
        /*
         * Wrapper for Atom-properties
         *
         * dereferences pointers and exposes reference-semantics
         */
        template<typename U, U* AtomViewBase::*p_prop, bool* AtomViewBase::*p_mod>
        class PropRef {
            // helper struct because std::bitset has no const_reference member type
            template<typename V>
            struct const_ref {
                using type = decltype(std::declval<const V>()[0]);
            };
        public:
            PropRef(AtomViewBase& at):at{at} {}
            // can only be explicitely constructed to reference an Atom
            PropRef(const PropRef&) = delete;
            // like a real reference, assigning changes the origin
            PropRef& operator=(const PropRef& rhs){
                *(at.*p_prop) = static_cast<const U&>(rhs);
                *(at.*p_mod) = true;
                return *this;
            }
            PropRef& operator=(const U& prop)
            {
                *(at.*p_prop) = prop;
                *(at.*p_mod) = true;
                return *this;
            }
            // convert to reference
            operator U&() {*(at.*p_mod) = true; return *(at.*p_prop);}
            operator const U&() const {return *(at.*p_prop);}
            // enable c_str-conversion for name
            template <typename R=U>
            const typename R::value_type* c_str() const {return (at.*p_prop)->c_str();}
            // enable array-access for Vec
            template <typename R=U>
            typename R::reference operator [](std::size_t i) {
                *(at.*p_mod) = true;
                return (*(at.*p_prop))[i];
            }
            template <typename R=U>
            typename const_ref<R>::type operator [](std::size_t i) const {
                return (*(at.*p_prop))[i];
            }
            // Comparison
            bool operator==(const U& rhs) const
            {
                return *(at.*p_prop) == rhs;
            }
            // Allow pointer-syntax for pse-pointer
            template <typename R=U, typename =typename std::enable_if<std::is_pointer<R>::value>::type>
            R operator->() const{
                return *(at.*p_prop);
            }
            // Allow member-acces via pointer-syntax for properties
            template <typename R=U, typename =typename std::enable_if<std::is_class<R>::value>::type>
            R* operator->(){
                *(at.*p_mod) = true;
                return at.*p_prop;
            }
            template <typename R=U, typename =typename std::enable_if<std::is_class<R>::value>::type>
            R const* operator->() const{
                return at.*p_prop;
            }
            // I/O
            friend std::ostream& operator<<(std::ostream& s, const PropRef& pr){
                s << static_cast<const U&>(pr);
                return s;
            }
            friend std::istream& operator>>(std::istream& s, PropRef& pr){
                s >> static_cast<U&>(pr);
                return s;
            }
        private:
            AtomViewBase& at;
        };
    };

    // Main interface declarations
    using Atom = AtomViewBase<void>;
    using constAtom = AtomViewBase<const void>;
}

#endif // ATOM_H
