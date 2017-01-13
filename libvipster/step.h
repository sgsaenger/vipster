#ifndef STEP_H
#define STEP_H

#include <config.h>
#include <vec.h>
#include <atom.h>
#include <bond.h>
#include <set>
#include <memory>

namespace Vipster{
class Step
{
public:
    Step();
    Step(const std::shared_ptr<PseMap> &pse);
    friend std::ostream& operator<< (std::ostream& s, const Step& st);
    void    newAtom(const Atom &at);                          //copy of atom
    void    newAtom(Atom&& at={});                            //move of atom
    void    newAtom(Atom at, AtomFmt fmt);                    //copy of atom (possibly too many)
    void    newAtoms(const std::vector<Atom> &v);             //batch copy
    void    delAtom(size_t idx);                              //delete
    void    setAtom(size_t idx,const Atom& at);               //replace with copy
    void    setAtom(size_t idx,Atom&& at);                    //replace with move
    void    setAtom(size_t idx,Atom at,AtomFmt fmt);          //replace with copy
    const Atom& getAtom(size_t idx)const;                     //get reference (const,bohr)
    Atom    getAtom(size_t idx, AtomFmt fmt);                 //get copy (formatted)
    const std::vector<Atom>& getAtoms(void) const noexcept;   //get const reference (bohr)
    std::vector<Atom> getAtoms(AtomFmt fmt);                  //get copy (formatted)
    size_t  getNat(void) const noexcept;                      //get number of atoms
    void    setCellDim(float cdm, bool scale=false, AtomFmt fmt=AtomFmt::Bohr);
    float   getCellDim(AtomFmt fmt=AtomFmt::Bohr) const noexcept;
    void    setCellVec(float v11, float v12, float v13,
                       float v21, float v22, float v23,
                       float v31, float v32, float v33,bool scale=false);
    void    setCellVec(const Vec &v1, const Vec &v2, const Vec &v3, bool scale=false);
    void    setCellVec(const Mat &vec, bool scale=false);
    Vec   getCenter(bool com=false) const noexcept;
    const Mat& getCellVec(void) const noexcept;
    std::set<std::string> getTypes(void)const noexcept;
    size_t  getNtyp(void) const noexcept;
    void  setComment(const std::string &s);
    const std::string& getComment(void) const noexcept;
    const std::vector<Bond>& getBonds() const;
    const std::vector<Bond>& getBonds(float cutfac) const;
    const std::vector<Bond>& getBondsCell() const;
    const std::vector<Bond>& getBondsCell(float cutfac) const;
    std::shared_ptr<PseMap> pse;
private:
    Atom formatAtom(Atom at, AtomFmt source, AtomFmt target);
    std::vector<Atom> formatAtoms(std::vector<Atom> atoms, AtomFmt source, AtomFmt target);
    void setBonds(float cutfac) const;
    void setBondsCell(float cutfac) const;
    enum class BondLevel { None, Molecule, Cell };
    //DATA following:
    std::string comment;
    std::vector<Atom> atoms;
    float celldim {1.};
    Mat cellvec {{ {{1.,0.,0.}}, {{0.,1.,0.}}, {{0.,0.,1.}} }};
    Mat invvec {{ {{1.,0.,0.}}, {{0.,1.,0.}}, {{0.,0.,1.}} }};
    mutable bool bonds_outdated{true};
    mutable BondLevel bonds_level{BondLevel::None};
    mutable float bondcut_factor{1.1};
    mutable std::vector<Bond> bonds;
};
std::ostream& operator<< (std::ostream& s, const Step& st);

}
#endif // STEP_H
