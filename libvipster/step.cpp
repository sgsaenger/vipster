#include "step.h"

using namespace Vipster;

Step::Step(AtomFmt fmt, const std::string &c)
    : StepMutable{std::make_shared<PeriodicTable>(),
                  std::make_shared<AtomList>(fmt),
                  std::make_shared<BondList>(),
                  std::make_shared<CellData>(),
                  std::make_shared<std::string>(c)}
{}

Step::Step(const Step& s)
    : StepMutable{s.pte,
                  std::make_shared<AtomList>(*s.atoms),
                  std::make_shared<BondList>(*s.bonds),
                  std::make_shared<CellData>(*s.cell),
                  std::make_shared<std::string>(*s.comment)}
{}

Step::Step(Step&& s)
    : StepMutable{s.pte,
                  s.atoms, s.bonds,
                  s.cell, s.comment}
{}

Step& Step::operator=(const Step& s)
{
    pte = s.pte;
    *atoms = *s.atoms;
    *bonds = *s.bonds;
    *cell = *s.cell;
    *comment = *s.comment;
    return *this;
}

Step& Step::operator=(Step&& s)
{
    pte = std::move(s.pte);
    atoms = std::move(s.atoms);
    bonds = std::move(s.bonds);
    cell = std::move(s.cell);
    comment = std::move(s.comment);
    return *this;
}

void Step::newAtom(std::string name, Vec coord, AtomProperties prop)
{
    AtomList& al = *atoms;
    al.coordinates.emplace_back(coord);
    al.elements.push_back(&*pte->find_or_fallback(name));
    al.properties.emplace_back(prop);
}

void Step::newAtoms(size_t i){
    size_t nat = getNat()+i;
    AtomList& al = *atoms;
    al.coordinates.resize(nat);
    al.elements.reserve(nat);
    for(size_t j=0; j<i; ++j){
        al.elements.push_back(&*pte->find_or_fallback(""));
    }
    al.properties.resize(nat);
}

void Step::delAtom(size_t _i){
    AtomList& al = *atoms;
    auto i = static_cast<long>(_i);
    al.coordinates.erase(al.coordinates.begin()+i);
    al.elements.erase(al.elements.begin()+i);
    al.properties.erase(al.properties.begin()+i);
}

void Step::enableCell(bool val) noexcept
{
    cell->enabled = val;
}

void Step::setCellVec(const Mat &vec, bool scale)
{
    Mat inv = Mat_inv(vec);
    cell->enabled = true;
    // 'scaling' means the system grows/shrinks with the cell
    if (scale == (atoms->fmt == AtomFmt::Crystal)){
        // scaled crystal or unscaled other coordinates stay as-is
        cell->cellvec = vec;
        cell->invvec = inv;
    }else{
        // unscaled crystal or scaled other coordinates have to be transformed
        auto tmpFmt = scale ? AtomFmt::Crystal : AtomFmt::Alat;
        auto tmp = formatAll(atoms->coordinates, atoms->fmt, tmpFmt);
        cell->cellvec = vec;
        cell->invvec = inv;
        atoms->coordinates = formatAll(tmp, tmpFmt, atoms->fmt);
    }
}

void Step::setCellDim(double cdm, CdmFmt fmt, bool scale)
{
    if(!(cdm>0)) {
        throw Error("Step::setCellDim(): "
                    "cell-dimension must be positive");
    }
    /*
     * 'scaling' means the system grows/shrinks with the cell
     * => relative coordinates stay the same
     */
    cell->enabled = true;
    bool relative = atomFmtRelative(atoms->fmt);
    if (scale != relative){
        double ratio;
        if (relative) {
            ratio = getCellDim(fmt) / cdm;
        } else {
            ratio = cdm / getCellDim(fmt);
        }
        for(auto& c: atoms->coordinates){
            c *= ratio;
        }
    }
    if (fmt == CdmFmt::Bohr) {
        cell->celldim = cdm * bohrrad;
    } else {
        cell->celldim = cdm;
    }
}

void Step::modWrap(){
    for(auto& at: asFmt(AtomFmt::Crystal)){
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
    auto oldNat = getNat();
    auto newNat = oldNat * fac;
    atoms->coordinates.reserve(newNat);
    atoms->elements.reserve(newNat);
    atoms->properties.reserve(newNat);
    auto handle = asFmt(AtomFmt::Crystal);
    auto cell = this->getCellVec();
    auto multiply = [&](uint8_t dir, size_t mult){
        cell[dir] *= mult;
        for(size_t i=1; i<mult; ++i){
            std::copy_n(atoms->coordinates.begin(), oldNat, std::back_inserter(atoms->coordinates));
            std::copy_n(atoms->elements.begin(), oldNat, std::back_inserter(atoms->elements));
            std::copy_n(atoms->properties.begin(), oldNat, std::back_inserter(atoms->properties));
            for(auto it = handle.begin()+i*oldNat; it!=handle.end(); ++it){
                it->coord[dir] += i;
            }
        }
    };
    if(x>1) multiply(0, x);
    if(y>1) multiply(1, y);
    if(z>1) multiply(2, z);
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
