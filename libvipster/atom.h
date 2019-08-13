#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <bitset>

#include "vec.h"
#include "periodictable.h"

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

    // supported formats | TODO: crystal/alat to the beginning so we can extend it more easily?
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
    template<bool isConst>
    class AtomViewBase{
        /* Wrapper classes that control access to the atoms
         *
         * should behave like references to the wrapped values,
         * but set the changed-flag whenever the values are accessed in a non-const way
         */
        struct _Vec{
            _Vec(AtomViewBase& at):at{at}{}
            // can only be explicitely constructed to reference an Atom
            _Vec(const AtomViewBase& at) =delete;
            // like real references, assigning changes the origin
            _Vec& operator=(const Vec& rhs){
                *at.coord_ptr = rhs;
                *at.coord_changed = true;
                return *this;
            }
            _Vec& operator=(const _Vec& rhs){
                *at.coord_ptr = static_cast<const Vec&>(rhs);
                *at.coord_changed = true;
                return *this;
            }
            // convert to reference
            operator Vec&() {*at.coord_changed = true; return *at.coord_ptr;}
            operator const Vec&() const {return *at.coord_ptr;}
            // array access
            Vec::value_type& operator[](std::size_t i) {
                *(at.coord_changed) = true;
                return (*at.coord_ptr)[i];
            }
            const Vec::value_type& operator[](std::size_t i) const {
                return (*at.coord_ptr)[i];
            }
            // comparison
            bool operator==(const Vec& rhs) const { return *at.coord_ptr == rhs;}
        private:
            AtomViewBase &at;
        };
        struct _Name{
            _Name(AtomViewBase& at):at{at}{}
            // can only be explicitely constructed to reference an Atom
            _Name(const AtomViewBase& at) =delete;
            // like real references, assigning changes the origin
            _Name& operator=(const std::string& rhs){
                *at.element_ptr = &*(*at.pte).find_or_fallback(rhs);
                *at.element_changed = true;
                return *this;
            }
            _Name& operator=(const _Name& rhs){
                *at.element_ptr = &*(*at.pte).find_or_fallback(rhs);
                *at.element_changed = true;
                return *this;
            }
            // convert to reference
            operator const std::string&() const {return (*at.element_ptr)->first;}
            // comparison
            bool operator==(const std::string& rhs) const { return (*at.element_ptr)->first == rhs;}
            // const char* access
            const char* c_str() const{
                return (*at.element_ptr)->first.c_str();
            }
            // I/O
            friend std::ostream& operator<<(std::ostream& s, const _Name &n) {
                return s << static_cast<const std::string&>(n);
            }
            friend std::istream& operator>>(std::istream& s, _Name &n){
                std::string tmp;
                s >> tmp;
                n = tmp;
                return s;
            }
        private:
            AtomViewBase &at;
        };
        struct _Element{
            _Element(AtomViewBase& at):at{at}{}
            // can only be explicitely constructed to reference an Atom
            _Element(const AtomViewBase& at) =delete;
            // Allow member-access via pointer-syntax
            const Element* operator->() const{return &(*at.element_ptr)->second;}
            // convert to reference
            operator const Element&() const {return (*at.element_ptr)->second;}
        private:
            AtomViewBase &at;
        };
        struct _Properties{
            _Properties(AtomViewBase& at):at{at}{}
            // can only be explicitely constructed to reference an Atom
            _Properties(const AtomViewBase& at) =delete;
            // like real references, assigning changes the origin
            _Properties& operator=(const AtomProperties& rhs){
                *at.prop_ptr = rhs;
                *at.prop_changed = true;
                return *this;
            }
            _Properties& operator=(const _Properties& rhs){
                *at.prop_ptr = static_cast<const AtomProperties&>(rhs);
                *at.prop_changed = true;
                return *this;
            }
            // convert to reference
            operator const AtomProperties&() const {return *at.prop_ptr;}
            // comparison
            bool operator==(const AtomProperties& rhs) const { return *at.prop_ptr == rhs; }
            // Allow pointer-syntax
            const AtomProperties* operator->() const{return at.prop_ptr;}
            AtomProperties* operator->() {*at.prop_changed = true; return at.prop_ptr;}
        private:
            AtomViewBase &at;
        };
    public:
        virtual ~AtomViewBase() = default;
        AtomViewBase() = delete;
        // copy constructor creates new object pointing to same data
        AtomViewBase(const AtomViewBase& at)
            : coord_ptr{at.coord_ptr},
              coord_changed{at.coord_changed},
              element_ptr{at.element_ptr},
              element_changed{at.element_changed},
              prop_ptr{at.prop_ptr},
              prop_changed{at.prop_changed},
              pte{at.pte}
        {}
        // allow Atom to constAtom conversion
        template<bool B> friend class AtomViewBase;
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomViewBase(const AtomViewBase<B>& at)
            : coord_ptr{at.coord_ptr},
              coord_changed{at.coord_changed},
              element_ptr{at.element_ptr},
              element_changed{at.element_changed},
              prop_ptr{at.prop_ptr},
              prop_changed{at.prop_changed},
              pte{at.pte}
        {}
        // copy assignment assigns values to old position
        AtomViewBase& operator=(const AtomViewBase& at){
            coord = at.coord;
            name = at.name;
            properties = at.properties;
            return *this;
        }

        // "Data members", encapsulated in wrapper objects
        _Vec        coord{*this};
        _Name       name{*this};
        _Element    type{*this};
        _Properties properties{*this};

        // comparison
        bool operator==(const AtomViewBase& rhs) const {
            return std::tie(name, coord, properties)
                   ==
                   std::tie(rhs.name, rhs.coord, rhs.properties);
        }
        bool operator!=(const AtomViewBase& rhs) const {
            return !operator==(rhs);
        }

    protected:
        /* actual constructors
         *
         * to be called from iterators, not standalone so far
         */
        AtomViewBase(Vec *co, bool *c_m, PeriodicTable::value_type **el, bool *e_m,
             AtomProperties *p, bool *p_m, PeriodicTable* pt)
            : coord_ptr{co}, coord_changed{c_m},
              element_ptr{el}, element_changed{e_m},
              prop_ptr{p}, prop_changed{p_m},
              pte{pt}
        {}
        // pointers to actual data-storage
        Vec                         *coord_ptr;
        bool                        *coord_changed;
        PeriodicTable::value_type*  *element_ptr;
        bool                        *element_changed;
        AtomProperties                  *prop_ptr;
        bool                        *prop_changed;
        /* pointer to source's periodic table
         *
         * using raw-pointer here because rest depends on life-time of parent, anyways
         */
        PeriodicTable               *pte;
    };

    // Main interface declarations
    using Atom = AtomViewBase<false>;
    using constAtom = AtomViewBase<true>;
}
#endif // ATOM_H
