#include "step.h"

using namespace Vipster;

Step::Step(AtomFmt fmt, const std::string &c)
    : StepMutable{std::make_shared<atom_source>(fmt),
                  std::make_shared<BondList>(),
                  std::make_shared<std::string>(c)}
{}

Step::Step(const Step& s)
    : StepMutable{std::make_shared<atom_source>(*s.atoms),
                  std::make_shared<BondList>(*s.bonds),
                  std::make_shared<std::string>(*s.comment)}
{
    for(auto& b: bonds->list){
        if(b.type){
            b.type = &*bonds->types.find(b.type->first);
        }
    }
}

Step::Step(Step&& s)
    : StepMutable{std::move(s)}
{}

Step& Step::operator=(const Step& s)
{
    *atoms = *s.atoms;
    *bonds = *s.bonds;
    for(auto& b: bonds->list){
        if(b.type){
            b.type = &*bonds->types.find(b.type->first);
        }
    }
    *comment = *s.comment;
    return *this;
}

Step& Step::operator=(Step&& s)
{
    atoms = std::move(s.atoms);
    bonds = std::move(s.bonds);
    comment = std::move(s.comment);
    return *this;
}

void Step::newAtom(std::string name, Vec coord, AtomProperties prop)
{
    atom_source& al = *atoms;
    al.coordinates.emplace_back(coord);
    al.elements.push_back(&*getPTE().find_or_fallback(name));
    al.properties.emplace_back(prop);
}

void Step::newAtoms(size_t i){
    size_t nat = getNat()+i;
    atom_source& al = *atoms;
    al.coordinates.resize(nat);
    al.elements.reserve(nat);
    for(size_t j=0; j<i; ++j){
        al.elements.push_back(&*getPTE().find_or_fallback(""));
    }
    al.properties.resize(nat);
}

void Step::delAtom(size_t _i){
    atom_source& al = *atoms;
    auto i = static_cast<long>(_i);
    al.coordinates.erase(al.coordinates.begin()+i);
    al.elements.erase(al.elements.begin()+i);
    al.properties.erase(al.properties.begin()+i);
}

void Step::setPTE(std::shared_ptr<PeriodicTable> newPTE)
{
    // reassign types
    for(auto &el: atoms->elements){
        el = &*newPTE->find_or_fallback(el->first);
    }
    // let context own new table
    atoms->ctxt.pte = std::move(newPTE);
}

void Step::setFmt(AtomFmt tgt, bool scale)
{
    if(tgt == this->atoms->ctxt.fmt){ return; }
    if(atomFmtRelative(tgt)){
        enableCell(true);
    }
    if(scale && (getNat() != 0)){
        auto tmp = asFmt(tgt);
        auto source = tmp.cbegin();
        auto target = begin();
        while(source != tmp.cend()){
            target->coord = source->coord;
            ++target;
            ++source;
        }
        atoms->ctxt.fmt = tgt;
    }
    atoms->ctxt.fmt = tgt;
}

void Step::enableCell(bool val) noexcept
{
    if(val && (atoms->ctxt.cell->matrix == Mat{})){
        // on initial enabling, fill-in matrix that fits current structure
        auto com = getCom(AtomFmt::Angstrom);
        setCellVec({{{{std::max(1., 2*com[0]+0.5), 0, 0}},
                     {{0, std::max(1., 2*com[1]+0.5), 0}},
                     {{0, 0, std::max(1., 2*com[2]+0.5)}}}});
    }else{
        atoms->ctxt.cell->enabled = val;
    }
}

void Step::setCellVec(const Mat &vec, bool scale)
{
    Mat inv = Mat_inv(vec);
    auto &cell = atoms->ctxt.cell;
    cell->enabled = true;
    bool isCrystal = atoms->ctxt.fmt == AtomFmt::Crystal;
    // 'scaling' means the system grows/shrinks with the cell
    if (scale == isCrystal){
        // scaled crystal or unscaled other coordinates stay as-is
        cell->matrix = vec;
        cell->inverse = inv;
    }else{
        // unscaled crystal or scaled other coordinates have to be transformed
        auto fmt = isCrystal ? cell->matrix * inv : cell->inverse * vec;
        for(auto &at: *this){
            at.coord = at.coord * fmt;
        }
        cell->matrix = vec;
        cell->inverse = inv;
    }
}

void Step::setCellDim(double cdm, AtomFmt fmt, bool scale)
{
    if(!atomFmtAbsolute(fmt)){
        throw Error{"Step::setCellDim: Invalid AtomFmt, needs to be absolute"};
    }
    if(!(cdm>0)) {
        throw Error("Step::setCellDim(): "
                    "cell-dimension must be positive");
    }
    /*
     * 'scaling' means the system grows/shrinks with the cell
     * => relative coordinates stay the same
     */
    enableCell(true);
    bool relative = atomFmtRelative(getFmt());
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
    atoms->ctxt.cell->dimension = cdm * detail::toAngstrom[static_cast<size_t>(fmt)];
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
    auto newNat = getNat() * fac;
    atoms->coordinates.reserve(newNat);
    atoms->elements.reserve(newNat);
    atoms->properties.reserve(newNat);
    auto cell = this->getCellVec();
    auto multiply = [&](uint8_t dir, size_t mult){
        auto oldNat = getNat();
        auto handle = asFmt(AtomFmt::Crystal);
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

void Step::modReshape(Mat newMat, double newDim, AtomFmt Fmt){
    auto oldCdm = getCellDim(Fmt);
    auto oldMat = getCellVec();
    if((newMat == oldMat) && (float_comp(newDim, oldCdm))){
        return;
    }
    modWrap();
    if(newMat == oldMat){
        // only changing cdm
        size_t fac = static_cast<size_t>(std::ceil(newDim/oldCdm));
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
        Vec newExtent = getExtent(newMat*newDim);
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
    setCellDim(newDim, Fmt);
    modCrop();
}
