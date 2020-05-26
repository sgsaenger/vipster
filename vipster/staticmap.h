#ifndef CUSTOMMAP_H
#define CUSTOMMAP_H

#include <map>

namespace Vipster {

template <typename Key, typename Value>
class StaticMap: protected std::map<Key, Value>
{
public:
    using map_t = std::map<Key, Value>;
    // access functions
    using map_t::begin;
    using map_t::end;
    using map_t::at;
    using map_t::find;
    using map_t::erase;
    using map_t::size;
    // types
    using typename map_t::iterator;
    using typename map_t::const_iterator;
    using typename map_t::key_type;
    using typename map_t::mapped_type;
    using typename map_t::value_type;
    using typename map_t::reference;
    // inherit constructors
    using map_t::map_t;
};

}

#endif // CUSTOMMAP_H
