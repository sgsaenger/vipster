#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <type_traits>
#include "vec.h"

//TODO: track changes in hidden (rename it, obviously)

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
            // like a real reference, constructing makes it point to the origin
            PropRef(const PropRef&) = default;
            // also like a real reference, assigning changes the origin
            PropRef& operator=(const PropRef& rhs){
                *p_prop = *(rhs.p_prop);
                *p_mod = true;
                return *this;
            }
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
             const FixVec *f, const uint8_t *h,
             const bool *c_m, const bool *p_m);
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
