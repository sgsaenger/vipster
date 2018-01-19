#ifndef ATOMPROPER_H
#define ATOMPROPER_H

#include "atom.h"

namespace Vipster {
    /*
     * Atom that owns its data
     */
    class AtomProper: public Atom{
    public:
        AtomProper(std::string name="", Vec coord={{0,0,0}}, float charge=0,
                   FixVec fix={{false,false,false}}, uint8_t hidden=false);
        AtomProper(const AtomProper& rhs);
        AtomProper& operator=(const AtomProper& rhs);
    private:
        std::string val_name;
        Vec val_coord;
        float val_charge;
        FixVec val_fix;
        uint8_t val_hidden;
        bool mod;
    };
}

#endif // ATOMPROPER_H
