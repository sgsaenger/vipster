#ifndef STEP_H
#define STEP_H

#include <global.h>
#include <config.h>
#include <vec.h>
#include <atom.h>
#include <bond.h>
#include <set>

namespace Vipster{

class Step
{
public:
    Step(PseMap pse);
    std::string comment;
    void    newAtom(std::string name="C",
                    Vec coord={0.,0.,0.},
                    float charge=0.,
                    std::array<bool,3> fix={false,false,false},
                    bool hidden=false,
                    Fmt fmt=Fmt::Bohr
    );                                                  //initialization of atom
    void    newAtom(const Atom &at);                  //copy of atom
    void    newAtom(Atom&& at);                       //move of atom
    void    newAtom(Atom at, Fmt fmt);                //copy of atom (possibly too many)
    void    newAtoms(size_t count);                     //batch creation
    void    delAtom(size_t idx);                        //delete
    void    setAtom(size_t idx,                         //modify/initialize
                    std::string name="C",
                    Vec coord={0.,0.,0.},
                    float charge=0.,
                    std::array<bool,3> fix={false,false,false},
                    bool hidden=false,
                    Fmt fmt=Fmt::Bohr
    );
    void    setAtom(size_t idx,const Atom& at);               //replace with copy
    void    setAtom(size_t idx,Atom&& at);                    //replace with move
    void    setAtom(size_t idx,Atom at,Fmt fmt);              //replace with copy
    const Atom& getAtom(size_t idx)const;                     //get reference (const,bohr)
    Atom    getAtomFmt(size_t idx, Fmt fmt);                  //get copy (formatted)
    const std::vector<Atom>& getAtoms(void) const noexcept;   //get const reference (bohr)
    std::vector<Atom> getAtomsFmt(Fmt fmt);                   //get copy (formatted)
    size_t  getNat(void) const noexcept;                      //get number of atoms
    void    setCellDim(float cdm, bool scale=false, Fmt fmt=Fmt::Bohr);
    float   getCellDim(Fmt fmt=Fmt::Bohr) const noexcept;
    void    setCellVec(float v11, float v12, float v13,
                       float v21, float v22, float v23,
                       float v31, float v32, float v33,bool scale=false);
    void    setCellVec(Vec v1, Vec v2, Vec v3,bool scale=false);
    void    setCellVec(std::array<Vec,3> vec,bool scale=false);
    Vec   getCenter(bool com=false) const;
    const std::array<Vec,3>& getCellVec(void) const noexcept;
    std::set<std::string> getTypes(void)const noexcept;
    size_t  getNtyp(void) const noexcept;
    const std::vector<Bond>& getBonds() const;
    const std::vector<Bond>& getBonds(float cutfac) const;
    const std::vector<Bond>& getBondsCell() const;
    const std::vector<Bond>& getBondsCell(float cutfac) const;
    mutable PseMap pse;
private:
    Atom formatAtom(Atom at, Fmt source, Fmt target);
    std::vector<Atom> formatAtoms(std::vector<Atom> atoms, Fmt source, Fmt target);
    void setBonds(float cutfac) const;
    void setBondsCell(float cutfac) const;
    enum class BondType { None, Molecule, Cell };
    //DATA following:
    std::vector<Atom> atoms;
    float celldim;
    std::array<Vec,3> cellvec;
    std::array<Vec,3> invvec;
    mutable bool bonds_outdated;
    mutable BondType bonds_level;
    mutable float bondcut_factor;
    mutable std::vector<Bond> bonds;
};


}
#endif // STEP_H
