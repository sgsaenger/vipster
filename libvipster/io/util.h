#ifndef IOUTIL_H
#define IOUTIL_H

#include "../io.h"

namespace Vipster::IO{

std::string trim(const std::string& str, const std::string& ws=" \t\r");
void intToCart(Step& s, const std::string& name, const std::array<size_t,3>&ids, Vec values);

}

#endif // IOUTIL_H
