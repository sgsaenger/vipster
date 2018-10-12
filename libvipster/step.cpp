#include "step.h"
#include "stepsel.h"

using namespace Vipster;

Step::Step(AtomFmt at_fmt, std::string comment)
    : StepMutable{std::make_shared<PseMap>(),
                  at_fmt,
                  std::make_shared<AtomList>(),
                  std::make_shared<BondList>(),
                  std::make_shared<CellData>(),
                  std::make_shared<std::string>(comment)}
{}

Step::Step(const Step& s)
    : StepMutable{s.pse, s.at_fmt,
                  std::make_shared<AtomList>(*s.atoms),
                  std::make_shared<BondList>(*s.bonds),
                  std::make_shared<CellData>(*s.cell),
                  std::make_shared<std::string>(*s.comment)}
{}

Step::Step(std::shared_ptr<PseMap> p, AtomFmt f,
           std::shared_ptr<AtomList> a, std::shared_ptr<BondList> b,
           std::shared_ptr<CellData> c, std::shared_ptr<std::string> s)
    : StepMutable{p,f,a,b,c,s}
{}

Step& Step::operator=(const Step& s)
{
    pse = s.pse;
    at_fmt = s.at_fmt;
    *atoms = *s.atoms;
    *bonds = *s.bonds;
    *cell = *s.cell;
    *comment = *s.comment;
    return *this;
}

Step Step::asFmt(AtomFmt tgt)
{
    return Step{pse, tgt, atoms, bonds, cell, comment};
}

StepSelection Step::select(std::string filter)
{
    return StepSelection{pse, at_fmt, std::make_shared<AtomSelection>(AtomSelection{{}, atoms}),
                         bonds, cell, comment, filter};
}

StepSelection Step::select(SelectionFilter filter)
{
    return StepSelection{pse, at_fmt, std::make_shared<AtomSelection>(AtomSelection{{}, atoms}),
                         bonds, cell, comment, filter};
}

StepSelConst Step::select(std::string filter) const
{
    return StepSelConst{pse, at_fmt, std::make_shared<AtomSelection>(AtomSelection{{}, atoms}),
                        bonds, cell, comment, filter};
}

StepSelConst Step::select(SelectionFilter filter) const
{
    return StepSelConst{pse, at_fmt, std::make_shared<AtomSelection>(AtomSelection{{}, atoms}),
                        bonds, cell, comment, filter};
}

void Step::newAtom(std::string name, Vec coord, AtomProperties prop)
{
    AtomList& al = *atoms;
    // Coordinates
    al.coordinates[static_cast<size_t>(at_fmt)].emplace_back(coord);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.names.emplace_back(name);
    al.name_changed = true;
    al.pse.push_back(&(*pse)[name]);
    // Properties
    al.properties.emplace_back(prop);
    al.prop_changed = true;
}

void Step::newAtom(const Atom& at){
    AtomList& al = *atoms;
    // Coordinates
    al.coordinates[static_cast<size_t>(at_fmt)].push_back(at.coord);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.names.push_back(at.name);
    al.name_changed = true;
    al.pse.push_back(&(*pse)[at.name]);
    // Properties
    al.properties.push_back(at.properties);
    al.prop_changed = true;
}

void Step::newAtoms(size_t i){
    size_t nat = getNat()+i;
    // Coordinates
    AtomList& al = *atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].resize(nat);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.names.resize(nat);
    al.name_changed = true;
    al.pse.resize(nat);
    // Properties
    al.properties.resize(nat);
    al.prop_changed = true;
}

void Step::delAtom(size_t _i){
    AtomList& al = *atoms;
    auto i = static_cast<long>(_i);
    // Coordinates
    al.coordinates[static_cast<size_t>(at_fmt)].erase(
        al.coordinates[static_cast<size_t>(at_fmt)].begin()+i);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.names.erase(al.names.begin()+i);
    al.name_changed = true;
    al.pse.erase(al.pse.begin()+i);
    // Properties
    al.properties.erase(al.properties.begin()+i);
    al.prop_changed = true;
}

void Step::setCellVec(const Mat &vec, bool scale)
{
    Mat inv = Mat_inv(vec);
    cell->enabled = true;
    evaluateCache();
    if (scale) {
        if (at_fmt == AtomFmt::Crystal) {
            /*
             * Crystal coordinates stay as-is
             */
            cell->cellvec = vec;
            cell->invvec = inv;
        } else {
            /*
             * All other cases need to be reformatted
             *
             * calculates crystal-coordinates as intermediate
             */
            asFmt(AtomFmt::Crystal).evaluateCache();
            cell->cellvec = vec;
            cell->invvec = inv;
            constexpr auto crystal = static_cast<size_t>(AtomFmt::Crystal);
            const auto target = static_cast<size_t>(at_fmt);
            atoms->coordinates[target] =
                    formatAll(atoms->coordinates[crystal],
                              AtomFmt::Crystal, at_fmt);
        }
    } else {
        if (at_fmt == AtomFmt::Crystal) {
            /*
             * Crystal needs to be reformatted
             *
             * Determine a valid buffer or use Alat
             */
            size_t buf = nAtFmt;
            for (size_t i=0; i<nAtFmt; ++i){
                if (i == static_cast<size_t>(AtomFmt::Crystal)){
                    continue;
                }
                if (!atoms->coord_outdated[i]){
                    buf = i;
                }
            }
            if (buf == nAtFmt) {
                buf = static_cast<size_t>(AtomFmt::Alat);
                asFmt(AtomFmt::Alat).evaluateCache();
            }
            cell->cellvec = vec;
            cell->invvec = inv;
            constexpr auto target = static_cast<size_t>(AtomFmt::Crystal);
            atoms->coordinates[target] =
                    formatAll(atoms->coordinates[buf],
                              static_cast<AtomFmt>(buf), AtomFmt::Crystal);
        } else {
            /*
             * All but crystal stay as-is
             */
            cell->cellvec = vec;
            cell->invvec = inv;
        }
    }
    atoms->coord_changed[static_cast<size_t>(at_fmt)] = true;
}

void Step::setCellDim(float cdm, CdmFmt fmt, bool scale)
{
    if(!(cdm>0)) {
        throw Error("Step::setCellDim(): "
                    "cell-dimension must be positive");
    }
    /*
     * 'scaling' means the systems grows/shrinks with the cell
     * => relative coordinates stay the same
     */
    cell->enabled = true;
    size_t int_fmt = static_cast<size_t>(at_fmt);
    bool relative = at_fmt>=AtomFmt::Crystal;
    if (scale != relative){
        float ratio;
        if (relative) {
            ratio = getCellDim(fmt) / cdm;
        } else {
            ratio = cdm / getCellDim(fmt);
        }
        for(auto& c: atoms->coordinates[int_fmt]){
            c *= ratio;
        }
    }
    atoms->coord_changed[int_fmt] = true;
    if (fmt == CdmFmt::Bohr) {
        cell->dimBohr = cdm;
        cell->dimAngstrom = cdm*bohrrad;
    } else {
        cell->dimAngstrom = cdm;
        cell->dimBohr = cdm*invbohr;
    }
}

template<>
size_t StepConst<AtomList>::getNat() const noexcept
{
    return atoms->names.size();
}

template<>
void StepConst<AtomList>::evaluateCache() const
{
    auto fmt = static_cast<size_t>(this->at_fmt);
    auto& atoms = *this->atoms;
    // if there are (only) local changes, invalidate other caches
    if(atoms.coord_changed[fmt]){
        for(size_t i=0; i<nAtFmt; ++i){
            if(i == fmt){
                atoms.coord_outdated[i] = false;
            }else if(atoms.coord_changed[i]){
                throw Error("Concurrent modification of a "
                            "single Step not supported.");
            }else{
                atoms.coord_outdated[i] = true;
            }
        }
        atoms.coord_changed[fmt] = false;
        bonds->outdated = true;
    } else {
        // if there are changes in (only one) other cache,
        // pull them in and invalidate the others
        size_t fmt_changed = nAtFmt;
        for(size_t i=0; i<nAtFmt; ++i){
            if (atoms.coord_changed[i]){
                if (fmt_changed == nAtFmt) {
                    fmt_changed = i;
                } else {
                    throw Error("Concurrent modification of a "
                                "single Step not supported.");
                }
            }
        }
        if (fmt_changed != nAtFmt) {
            atoms.coordinates[fmt] = formatAll(atoms.coordinates[fmt_changed],
                                               static_cast<AtomFmt>(fmt_changed),
                                               at_fmt);
            for (size_t i=0; i<nAtFmt; ++i){
                if((i == fmt) || (i == fmt_changed)){
                    atoms.coord_changed[i] = false;
                    atoms.coord_outdated[i] = false;
                }else{
                    atoms.coord_outdated[i] = true;
                }
            }
            bonds->outdated = true;
        } else if(atoms.coord_outdated[fmt]){
            // validate this cache if needed
            size_t fmt_uptodate = nAtFmt;
            for(size_t i=0; i<nAtFmt; ++i){
                if(!atoms.coord_outdated[i]){
                    fmt_uptodate = i;
                    break;
                }
            }
            if(fmt_uptodate == nAtFmt){
                throw Error("No valid atom-cache");
            }
            atoms.coordinates[fmt] = formatAll(atoms.coordinates[fmt_uptodate],
                                               static_cast<AtomFmt>(fmt_uptodate),
                                               at_fmt);
            atoms.coord_outdated[fmt] = false;
        }
    }
    // if some atom-types changed, update pse pointers
    if(atoms.name_changed){
        if(size_t nat = atoms.names.size()){
            for(size_t i=0; i<nat; ++i){
                atoms.pse[i] = &(*pse)[atoms.names[i]];
            }
        }
        atoms.name_changed = false;
        bonds->outdated = true;
    }
}
