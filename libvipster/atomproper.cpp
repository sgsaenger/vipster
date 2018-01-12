#include "atomproper.h"

using namespace Vipster;

AtomProper::AtomProper(std::string name, Vec coord, float charge,
                       FixVec fix, char hidden)
    : Atom{&val_name, &val_coord, &val_charge, &val_fix, &val_hidden, &mod},
      val_name{name}, val_coord{coord}, val_charge{charge},
      val_fix{fix}, val_hidden{hidden} {}

AtomProper::AtomProper(const AtomProper& rhs)
    : Atom{&val_name, &val_coord, &val_charge, &val_fix, &val_hidden, &mod},
      val_name{rhs.name}, val_coord(rhs.coord), val_charge{rhs.charge},
      val_fix(rhs.fix), val_hidden{rhs.hidden} {}

AtomProper& AtomProper::operator=(const AtomProper& rhs){
    val_name = rhs.val_name;
    val_coord = rhs.val_coord;
    val_charge = rhs.val_charge;
    val_fix = rhs.val_fix;
    val_hidden = rhs.val_hidden;
    return *this;
}
