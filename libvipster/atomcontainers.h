#ifndef ATOMCONTAINERS_H
#define ATOMCONTAINERS_H

#include "atom.h"
#include "vec.h"
#include <array>

namespace Vipster{

/*
 * Basic serial Atom container
 *
 * Stores atom in separate vectors
 */
struct AtomList{
    // Coordinates
    // one buffer per Vipster::AtomFmt
    std::array<std::vector<Vec>, nAtFmt> coordinates;
    std::array<bool, nAtFmt>             coord_changed;
    std::array<bool, nAtFmt>             coord_outdated;
    // Names (synced with type-pointers)
    std::vector<std::string>        names;
    bool                            name_changed;
    std::vector<PseEntry*>          pse;
    // Properties
    std::vector<AtomProperties>     properties;
    bool                            prop_changed;
};

/*
 * Selection container
 *
 * contains indices of selected atoms in AtomList
 */
struct AtomSelection{
    std::vector<size_t> indices;
    std::shared_ptr<AtomList> atoms;
};

/*
 * Iterator for serial Atom container
 */
template<typename T>
class AtomListIterator: private T
{
public:
    AtomListIterator(const std::shared_ptr<AtomList> &atoms,
                     AtomFmt fmt, size_t idx)
        : T{&atoms->coordinates[static_cast<uint8_t>(fmt)][idx],
            &atoms->coord_changed[static_cast<uint8_t>(fmt)],
            &atoms->names[idx],
            &atoms->name_changed,
            &atoms->properties[idx],
            &atoms->pse[idx],
            &atoms->prop_changed,
        }, atoms{atoms}, fmt{fmt}, idx{idx}
    {}
    AtomListIterator& operator++(){
        ++idx;
        ++(this->coord_ptr);
        ++(this->name_ptr);
        ++(this->prop_ptr);
        ++(this->pse_ptr);
        return *this;
    }
    AtomListIterator& operator+=(size_t i){
        idx += i;
        this->coord_ptr += i;
        this->name_ptr += i;
        this->prop_ptr += i;
        this->pse_ptr += i;
        return *this;
    }
    AtomListIterator operator+(size_t i){
        AtomListIterator copy{*this};
        return copy+=i;
    }
    T&      operator*() const {
        //const-ness of iterator is separate of const-ness of atoms!
        return static_cast<T&>(*const_cast<AtomListIterator*>(this));
    }
    T*      operator->() const {
        return &(operator*());
    }
    bool    operator==(const AtomListIterator& rhs) const noexcept{
        return (atoms == rhs.atoms) && (fmt == rhs.fmt) && (idx == rhs.idx);
    }
    bool    operator!=(const AtomListIterator& rhs) const noexcept{
        return !(*this == rhs);
    }
    size_t getIdx() const noexcept{
        return idx;
    }
private:
    std::shared_ptr<AtomList> atoms;
    AtomFmt fmt;
    size_t idx;
};


/*
 * Iterator for Atom selection
 *
 * dereferences selection-indices
 */
template<typename T>
class AtomSelIterator: private T
{
public:
    AtomSelIterator(const std::shared_ptr<AtomSelection> &selection,
                    AtomFmt fmt, size_t idx)
    //TODO: introduce a terminal-object (when c++17 is ready?)
        : T{&selection->atoms->coordinates[static_cast<size_t>(fmt)][idx<selection->indices.size()?selection->indices[idx]:0],
            &selection->atoms->coord_changed[static_cast<size_t>(fmt)],
            &selection->atoms->names[idx<selection->indices.size()?selection->indices[idx]:0],
            &selection->atoms->name_changed,
            &selection->atoms->properties[idx<selection->indices.size()?selection->indices[idx]:0],
            &selection->atoms->pse[idx<selection->indices.size()?selection->indices[idx]:0],
            &selection->atoms->prop_changed},
          selection{selection}, fmt{fmt}, idx{idx}
    {}
    AtomSelIterator& operator++(){
        ++idx;
        auto diff = selection->indices[idx] - selection->indices[idx-1];
        this->coord_ptr += diff;
        this->name_ptr += diff;
        this->prop_ptr += diff;
        this->pse_ptr += diff;
        return *this;
    }
    AtomSelIterator& operator+=(size_t i){
        idx += i;
        auto diff = selection->indices[idx] - selection->indices[idx-i];
        this->coord_ptr += diff;
        this->name_ptr += diff;
        this->prop_ptr += diff;
        this->pse_ptr += diff;
        return *this;
    }
    AtomSelIterator operator+(size_t i){
        AtomSelIterator copy{*this};
        return copy+=i;
    }
    T&  operator*() const {
        return static_cast<T&>(*const_cast<AtomSelIterator*>(this));
    }
    T*  operator->() const {
        return &(operator*());
    }
    bool    operator==(const AtomSelIterator& rhs) const noexcept{
        return (selection == rhs.selection) && (fmt == rhs.fmt) && (idx == rhs.idx);
    }
    bool    operator!=(const AtomSelIterator& rhs) const noexcept{
        return !(*this == rhs);
    }
    size_t getIdx() const noexcept{
        return idx;
    }
private:
    std::shared_ptr<AtomSelection> selection;
    AtomFmt fmt;
    size_t idx;
};

}

#endif // ATOMCONTAINERS_H
