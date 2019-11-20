#include "step.h"
#include "stepsel.h"

using namespace Vipster;

Step::Step(AtomFmt fmt, const std::string &c)
    : StepMutable{std::make_shared<PeriodicTable>(),
                  fmt,
                  std::make_shared<AtomList>(),
                  std::make_shared<BondList>(),
                  std::make_shared<CellData>(),
                  std::make_shared<std::string>(c)}
{}

Step::Step(const Step& s)
    : StepMutable{s.pte, s.at_fmt,
                  std::make_shared<AtomList>(*s.atoms),
                  std::make_shared<BondList>(*s.bonds),
                  std::make_shared<CellData>(*s.cell),
                  std::make_shared<std::string>(*s.comment)}
{}

Step::Step(Step&& s)
    : StepMutable{s.pte, s.at_fmt,
                  s.atoms, s.bonds,
                  s.cell, s.comment}
{}

Step::Step(std::shared_ptr<PeriodicTable> p, AtomFmt f,
           std::shared_ptr<AtomList> a, std::shared_ptr<BondList> b,
           std::shared_ptr<CellData> c, std::shared_ptr<std::string> s)
    : StepMutable{p,f,a,b,c,s}
{}

Step& Step::operator=(const Step& s)
{
    pte = s.pte;
    at_fmt = s.at_fmt;
    *atoms = *s.atoms;
    *bonds = *s.bonds;
    *cell = *s.cell;
    *comment = *s.comment;
    return *this;
}

Step& Step::operator=(Step&& s)
{
    pte = std::move(s.pte);
    at_fmt = s.at_fmt;
    atoms = std::move(s.atoms);
    bonds = std::move(s.bonds);
    cell = std::move(s.cell);
    comment = std::move(s.comment);
    return *this;
}

Step Step::asFmt(AtomFmt tgt)
{
    auto tmp = Step{pte, tgt, atoms, bonds, cell, comment};
    tmp.evaluateCache();
    return tmp;
}

void Step::newAtom(std::string name, Vec coord, AtomProperties prop)
{
    AtomList& al = *atoms;
    // Coordinates
    al.coordinates[static_cast<size_t>(at_fmt)].emplace_back(coord);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.elements.push_back(&*pte->find_or_fallback(name));
    al.element_changed = true;
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
    al.elements.push_back(&*pte->find_or_fallback(at.name));
    al.element_changed = true;
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
    al.elements.reserve(nat);
    for(size_t j=0; j<i; ++j){
        al.elements.push_back(&*pte->find_or_fallback(""));
    }
    al.element_changed = true;
    // Properties
    al.properties.resize(nat);
    al.prop_changed = true;
}

void Step::newAtoms(const AtomList& at){
    size_t nat = getNat() + at.elements.size();
    // Coordinates
    AtomList& al = *this->atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].reserve(nat);
    al.coordinates[static_cast<size_t>(at_fmt)].insert(
                al.coordinates[static_cast<size_t>(at_fmt)].end(),
                at.coordinates[static_cast<size_t>(at_fmt)].begin(),
                at.coordinates[static_cast<size_t>(at_fmt)].end());
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.elements.reserve(nat);
    for(const auto &el: at.elements){
        al.elements.push_back(&*pte->find_or_fallback(el->first));
    }
    al.element_changed = true;
    // Properties
    al.properties.reserve(nat);
    al.properties.insert(al.properties.end(), at.properties.begin(), at.properties.end());
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
    al.elements.erase(al.elements.begin()+i);
    al.element_changed = true;
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
            // evaluate crystal-cache via side-effect
            asFmt(AtomFmt::Crystal);
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
                // evaluate alat-cache via side-effect;
                asFmt(AtomFmt::Alat);
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

void Step::setCellDim(double cdm, CdmFmt fmt, bool scale)
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
        double ratio;
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

size_t AtomList::getNat() const noexcept
{
    return elements.size();
}

void AtomList::evaluateCache(const StepConst<AtomList> &step)
{
    auto fmt = static_cast<size_t>(step.at_fmt);
    // if there are (only) local changes, invalidate other caches
    if(coord_changed[fmt]){
        for(size_t i=0; i<nAtFmt; ++i){
            if(i == fmt){
                coord_outdated[i] = false;
            }else if(coord_changed[i]){
                throw Error("Concurrent modification of a "
                            "single Step not supported.");
            }else{
                coord_outdated[i] = true;
            }
        }
        coord_changed[fmt] = false;
        step.bonds->outdated = true;
    } else {
        // if there are changes in (only one) other cache,
        // pull them in and invalidate the others
        size_t fmt_changed = nAtFmt;
        for(size_t i=0; i<nAtFmt; ++i){
            if (coord_changed[i]){
                if (fmt_changed == nAtFmt) {
                    fmt_changed = i;
                } else {
                    throw Error("Concurrent modification of a "
                                "single Step not supported.");
                }
            }
        }
        if (fmt_changed != nAtFmt) {
            coordinates[fmt] = step.formatAll(coordinates[fmt_changed],
                                              static_cast<AtomFmt>(fmt_changed),
                                              step.at_fmt);
            for (size_t i=0; i<nAtFmt; ++i){
                if((i == fmt) || (i == fmt_changed)){
                    coord_changed[i] = false;
                    coord_outdated[i] = false;
                }else{
                    coord_outdated[i] = true;
                }
            }
            step.bonds->outdated = true;
        } else if(coord_outdated[fmt]){
            // validate this cache if needed
            size_t fmt_uptodate = nAtFmt;
            for(size_t i=0; i<nAtFmt; ++i){
                if(!coord_outdated[i]){
                    fmt_uptodate = i;
                    break;
                }
            }
            if(fmt_uptodate == nAtFmt){
                throw Error("No valid atom-cache");
            }
            coordinates[fmt] = step.formatAll(coordinates[fmt_uptodate],
                                              static_cast<AtomFmt>(fmt_uptodate),
                                              step.at_fmt);
            coord_outdated[fmt] = false;
        }
    }
    // if some atom-types changed, bonds may be invalid
    if(element_changed){
        element_changed = false;
        step.bonds->outdated = true;
    }
}

void Step::modWrap(){
    for(Atom& at: asFmt(AtomFmt::Crystal)){
        at.coord[0] -= std::floor(at.coord[0]);
        at.coord[1] -= std::floor(at.coord[1]);
        at.coord[2] -= std::floor(at.coord[2]);
    }
}

void Step::modCrop(){
    std::vector<size_t> toRemove;
    toRemove.reserve(this->getNat());
    const auto handle = asFmt(AtomFmt::Crystal);
    for(auto it=handle.cbegin(); it!=handle.cend(); ++it){
        if((it->coord[0]>=1) | (it->coord[0]<0) |
           (it->coord[1]>=1) | (it->coord[1]<0) |
           (it->coord[2]>=1) | (it->coord[2]<0)){
            toRemove.push_back(it.getIdx());
        }
    }
    for(auto it=toRemove.rbegin(); it!=toRemove.rend(); ++it){
        delAtom(*it);
    }
}

void Step::modMultiply(size_t x, size_t y, size_t z){
    auto fac = x*y*z;
    if(fac == 0){
        throw Error("Cannot eradicate atoms via modMultiply");
    }else if(fac == 1){
        return;
    }
    auto handle = asFmt(AtomFmt::Crystal);
    auto cell = this->getCellVec();
    auto multiply = [&](uint8_t dir, size_t mult){
        auto atoms = handle.getAtoms();
        auto oldNat = handle.getNat();
        cell[dir] *= mult;
        for(size_t i=1; i<mult; ++i){
            handle.newAtoms(atoms);
            auto refIt = handle.begin();
            for(auto it=refIt+i*oldNat; it!=refIt+(i+1)*oldNat; ++it){
                it->coord[dir] += i;
            }
        }
    };
    if(x>1){
        multiply(0, x);
    }
    if(y>1){
        multiply(1, y);
    }
    if(z>1){
        multiply(2, z);
    }
    setCellVec(cell);
}

void Step::modAlign(uint8_t step_dir, uint8_t target_dir){
    auto target = Vec{};
    target.at(target_dir) = 1;
    auto source = getCellVec().at(step_dir);
    source /= Vec_length(source);
    if(target == source){
        return;
    }
    auto axis = Vec_cross(source, target);
    axis /= Vec_length(axis);
    auto cos = Vec_dot(source, target);
    auto icos = 1-cos;
    auto sin = -std::sqrt(1-cos*cos);
    Mat rotMat = {Vec{icos * axis[0] * axis[0] + cos,
                      icos * axis[0] * axis[1] - sin * axis[2],
                      icos * axis[0] * axis[2] + sin * axis[1]},
                  Vec{icos * axis[1] * axis[0] + sin * axis[2],
                      icos * axis[1] * axis[1] + cos,
                      icos * axis[1] * axis[2] - sin * axis[0]},
                  Vec{icos * axis[2] * axis[0] - sin * axis[1],
                      icos * axis[2] * axis[1] + sin * axis[0],
                      icos * axis[2] * axis[2] + cos}};
    Mat oldCell = this->getCellVec();
    Mat newCell = oldCell*rotMat;
    setCellVec(newCell, true);
}

void Step::modReshape(Mat newMat, double newCdm, CdmFmt cdmFmt){
    auto oldCdm = getCellDim(cdmFmt);
    auto oldMat = getCellVec();
    if((newMat == oldMat) && (float_comp(newCdm, oldCdm))){
        return;
    }
    modWrap();
    size_t fac;
    if(newMat == oldMat){
        // only changing cdm
        fac = static_cast<size_t>(std::ceil(newCdm/oldCdm));
        modMultiply(fac, fac, fac);
    }else{
        // change vectors or both
        auto getMin = [](const Mat& m){
            return Vec{std::min({m[0][0], m[1][0], m[2][0]}),
                       std::min({m[0][1], m[1][1], m[2][1]}),
                       std::min({m[0][2], m[1][2], m[2][2]})};
        };
        auto oldMin = getMin(oldMat);
        auto newMin = getMin(newMat);
        auto getExtent = [](const Mat& m){
            return Vec{std::abs(m[0][0]) + std::abs(m[1][0]) + std::abs(m[2][0]),
                       std::abs(m[0][1]) + std::abs(m[1][1]) + std::abs(m[2][1]),
                       std::abs(m[0][2]) + std::abs(m[1][2]) + std::abs(m[2][2])
                       };
        };
        Vec newExtent = getExtent(newMat*newCdm);
        Vec oldExtent = getExtent(oldMat*oldCdm);
        size_t fac[3] = {static_cast<size_t>(std::ceil(newExtent[0]/oldExtent[0])),
                         static_cast<size_t>(std::ceil(newExtent[1]/oldExtent[1])),
                         static_cast<size_t>(std::ceil(newExtent[2]/oldExtent[2]))};
        bool off[3]{};
        for(size_t i=0; i<3; ++i){
            if(newMin[i]<oldMin[i]){
                fac[i] += 1;
                off[i] = true;
            }
        }
        modMultiply(fac[0], fac[1], fac[2]);
        if(off[0] || off[1] || off[2]){
            Vec offset{
                off[0] ? -1./fac[0] : 0.,
                off[1] ? -1./fac[1] : 0.,
                off[2] ? -1./fac[2] : 0.,
            };
            modShift(offset);
        }
    }
    setCellVec(newMat);
    setCellDim(newCdm, cdmFmt);
    modCrop();
}
