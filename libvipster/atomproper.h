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
                   FixVec fix={{false,false,false}}, char hidden=false);
        AtomProper(const AtomProper& rhs);
        AtomProper& operator=(const AtomProper& rhs);
    private:
        std::string val_name;
        Vec val_coord;
        float val_charge;
        FixVec val_fix;
        char val_hidden;
        bool mod;
    };
}

#endif // ATOMPROPER_H
