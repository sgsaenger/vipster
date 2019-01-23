#include "util.h"

using namespace Vipster;

std::string IO::trim(const std::string& str, const std::string& ws){
    const auto strBegin = str.find_first_not_of(ws);
    if(strBegin == std::string::npos){
        return "";
    }
    const auto strEnd = str.find_last_not_of(ws);
    return str.substr(strBegin, strEnd-strBegin+1);
}

void IO::intToCart(Step& s, const std::string& name, const std::array<size_t,3>&ids, Vec values)
{
    if(!ids[0]){
        // cartesian
        s.newAtom(name, values);
        return;
    }
    if(ids[1]){
        // angle
        values[1] *= deg2rad;
    }else{
        // shift in x-direction
        s.newAtom(name, s[ids[0]-1].coord + Vec{values[0],0,0});
        return;
    }
    if(ids[2]){
        // dihedral
        values[2] *= deg2rad;
    }else{
        // rotate in xz-plane TODO IS THAT CORRECT PLZ CHECK ME KTHXBYE
        s.newAtom(name, s[ids[0]-1].coord + Vec{
            std::cos(values[1])*values[0],
            0,
            std::sin(values[1])*values[0]});
        return;
    }
    // full dihedral-containing rotation
    auto v1 = s[ids[0]-1].coord - s[ids[1]-1].coord;
    auto v2 = s[ids[0]-1].coord - s[ids[2]-1].coord;
    auto n1 = Vec_cross(v1, v2);
    auto n2 = Vec_cross(v1, n1);
    n1 = n1 * -std::sin(values[2]) / Vec_length(n1);
    n2 = n2 * std::cos(values[2]) / Vec_length(n2);
    auto v3 = n1 + n2;
    v3 = v3 * values[0] * std::sin(values[1]) / Vec_length(v3);
    v1 = v1 * values[0] * std::cos(values[1]) / Vec_length(v1);
    s.newAtom(name, s[ids[0]-1].coord + v3 - v1);
}
