#ifndef NAMEDENUM_H
#define NAMEDENUM_H

#include "staticmap.h"
#include <vector>
#include <string>

namespace Vipster{
class NamedEnum: public StaticMap<int, std::string>
{
public:
    NamedEnum(int value, const std::vector<std::string> &names);
    operator int() const;
    operator const std::string&() const;
    int value() const;
    const std::string& name() const;
    std::vector<std::string> names() const;
    NamedEnum& operator=(int);
    NamedEnum& operator=(const std::string&);
private:
    int val;
};

}

#endif // NAMEDENUM_H
