#ifndef IOUTIL_H
#define IOUTIL_H

#include "fileio.h"

namespace Vipster{

std::string trim(const std::string& str, const std::string& ws=" \t\r");
std::pair<std::string, std::string> stripComment(const std::string& str, const std::string& m="#");
void intToCart(Step& s, const std::string& name, const std::array<size_t,3>&ids, Vec values);

}

#endif // IOUTIL_H
