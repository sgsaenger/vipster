#ifndef STEP_H
#define STEP_H

#include "config.h"
#include "definitions.h"

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
    Atom  getAtomFmt(size_t idx, Fmt fmt);                    //get copy (formatted)
    const std::vector<Atom>& getAtoms(void) const noexcept;   //get const reference (bohr)
    std::vector<Atom> getAtomsFmt(Fmt fmt);                   //get copy (formatted)
    size_t  getNat(void) const noexcept;                        //get number of atoms
    void    setCellDim(float cdm, bool scale=false);
    float   getCellDim(void) const noexcept;
    void    setCellVec(float v11, float v12, float v13,
                       float v21, float v22, float v23,
                       float v31, float v32, float v33,bool scale=false);
    void    setCellVec(Vec v1, Vec v2, Vec v3,bool scale=false);
    void    setCellVec(std::array<Vec,3> vec,bool scale=false);
    Vec   getCenter(bool com=false) const;
    const std::array<Vec,3>& getCellVec(void) const noexcept;
    std::set<std::string> getTypes(void)const noexcept;
    size_t  getNtyp(void) const noexcept;
    const std::array<std::vector<Bond>,8>& getBonds() const;
    const std::array<std::vector<Bond>,8>& getBonds(float cutfac) const;
    std::array<std::vector<std::array<Vec,2>>,8> getBondOffsets() const;
    mutable PseMap pse;
private:
    Atom formatAtom(Atom at, Fmt source, Fmt target);
    std::vector<Atom> formatAtoms(std::vector<Atom> atoms, Fmt source, Fmt target);
    void setBonds(float cutfac) const;
    //DATA following:
    std::vector<Atom> atoms;
    float celldim;
    std::array<Vec,3> cellvec;
    std::array<Vec,3> invvec;
    mutable bool bonds_outdated;
    mutable float bondcut_factor;
    mutable std::array<std::vector<Bond>,8> bonds;
};


}
#endif // STEP_H
