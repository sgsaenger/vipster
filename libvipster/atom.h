#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <type_traits>
#include "vec.h"

//TODO: track changes in hidden (rename it, obviously)

/*
 * use const where possible!
 * non-const-access may trigger reevaluation of step-properties depending
 * on atom-properties (Bonds!)
 *
 * current architecture may prove problematic when only few atoms in step
 * need to be modified
 * - non-const reading triggers unnecessary mod-flag-setting
 * - but needs to be non-const for real modifications
 * workaround: const-iteration with non-const Step::operator[] access?
 * needs user-awareness -> bad!
 *
 * TODO TODO TODO
 *
 * but will stay like this for now...
 */

/*
 * TODO:
 *
 * Benchmark performance, can references reduce overhead?
 * Or virtual functions after all?
 * References prohibit copy-construction (or just assignment?), problem?
 * Move-assign/construct PropRef should be possible
 */

namespace Vipster{
    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };

    using FixVec = std::array<bool, 3>;

    /*
     * Basic atom interface
     *
     * Provides access-wrappers to properties
     */
    class Atom{
    public:
        template<typename T, bool Name = false>
        class PropRef {
            friend class Atom;
        public:
            PropRef(const T *prop, const bool *mod)
                : p_prop{const_cast<T*>(prop)},
                  p_mod{const_cast<bool*>(mod)} {}
            PropRef(const PropRef&) = default;
            PropRef& operator=(const PropRef& rhs){
                *p_prop = *(rhs.p_prop);
                *p_mod = true;
                return *this;
            }
            // TODO! EXPENSIVE! use const Atom for read-only access!
            operator T&() {*p_mod = true; return *p_prop; }
            operator const T&() const { return *p_prop; }
            PropRef<T>& operator=(const T& prop)
            {
                *p_prop = prop;
                *p_mod = true;
                return *this;
            }
            template <typename R=T, typename = std::enable_if<Name>>
            const typename R::value_type* c_str() const {return p_prop->c_str();}
            // TODO! EXPENSIVE! use const Atom for read-only access!
            template <typename R=T, typename = std::enable_if<std::is_array<T>::value>>
            typename R::reference operator [](size_t i) {*p_mod = true; return (*p_prop)[i];}
            template <typename R=T, typename = std::enable_if<std::is_array<T>::value>>
            typename R::const_reference operator [](size_t i) const {return (*p_prop)[i];}
            bool operator==(const PropRef<T>&rhs) const noexcept
            {
                return *p_prop == *(rhs.p_prop);
            }
        private:
            T* p_prop;
            bool* p_mod;
        };
        virtual ~Atom() = default;
        PropRef<std::string> name;
        PropRef<Vec> coord;
        PropRef<float> charge;
        PropRef<FixVec> fix;
        PropRef<uint8_t> hidden;
    protected:
        Atom(const std::string *n, const Vec *co, const float *ch,
             const FixVec *f, const uint8_t *h, const bool *m);
        Atom& operator++();
    };
    //string may fail to be converted implicitely
    std::ostream& operator<<(std::ostream& s, const Atom::PropRef<std::string>& pr);
    std::istream& operator>>(std::istream& s, Atom::PropRef<std::string>& pr);
    bool operator==(const std::string& s, const Atom::PropRef<std::string>& pr);
    bool operator==(const Atom::PropRef<std::string>& pr, const std::string& s);

    bool operator==(const Atom &a1, const Atom &a2);
    bool operator!=(const Atom &a1, const Atom &a2);
}

#endif // ATOM_H
