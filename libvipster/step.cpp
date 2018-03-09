#include "step.h"

using namespace Vipster;

void Step::enableCell(bool b) noexcept
{
    cell->enabled = b;
}

size_t Step::getNat() const noexcept{
    evaluateCache();
    return atoms->names.size();
}

void Step::newAtom(){
    evaluateCache();
    AtomList& al = *atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].emplace_back();
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    al.names.emplace_back();
    al.charges.emplace_back();
    al.properties.emplace_back();
    al.pse.push_back(&(*pse)[""]);
}

void Step::newAtom(std::string name, Vec coord, float charge,
                   std::bitset<nAtProp> prop)
{
    evaluateCache();
    AtomList& al = *atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].emplace_back(coord);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    al.names.emplace_back(name);
    al.charges.emplace_back(charge);
    al.properties.emplace_back(prop);
    al.pse.push_back(&(*pse)[name]);
}

void Step::newAtom(const Atom& at){
    evaluateCache();
    AtomList& al = *atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].push_back(at.coord);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    al.names.push_back(at.name);
    al.charges.push_back(at.charge);
    al.properties.push_back(at.properties);
    al.pse.push_back(&(*pse)[at.name]);
}

void Step::newAtoms(size_t i){
    evaluateCache();
    size_t nat = getNat()+i;
    AtomList& al = *atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].resize(nat);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    al.names.resize(nat);
    al.charges.resize(nat);
    al.properties.resize(nat);
    al.pse.resize(nat);
}

void Step::delAtom(long i){
    evaluateCache();
    AtomList& al = *atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].erase(
        al.coordinates[static_cast<size_t>(at_fmt)].begin()+i);
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    al.names.erase(al.names.begin()+i);
    al.charges.erase(al.charges.begin()+i);
    al.properties.erase(al.properties.begin()+i);
    al.pse.erase(al.pse.begin()+i);
}

Step::iterator Step::begin() noexcept {
    evaluateCache();
    return {atoms, at_fmt, 0};
}

Step::constIterator Step::begin() const noexcept {
    evaluateCache();
    return {atoms, at_fmt, 0};
}

Step::constIterator Step::cbegin() const noexcept {
    evaluateCache();
    return {atoms, at_fmt, 0};
}

Step::iterator Step::end() noexcept {
    return {atoms, at_fmt, getNat()};
}

Step::constIterator Step::end() const noexcept {
    return {atoms, at_fmt, getNat()};
}

Step::constIterator Step::cend() const noexcept {
    return {atoms, at_fmt, getNat()};
}

Atom Step::operator[](size_t i) {
    evaluateCache();
    return *iterator{atoms, at_fmt, i};
}

void Step::setCellVec(const Mat &vec, bool scale)
{
    Mat inv = Mat_inv(vec);
    enableCell(true);
    if (scale) {
        if (at_fmt == AtomFmt::Crystal) {
            /*
             * Crystal coordinates stay as-is
             * but bonds should be reevaluated
             */
            cell->cellvec = vec;
            cell->invvec = inv;
            bonds->outdated = true;
        } else {
            /*
             * All other cases need to be reformatted
             *
             * calculates crystal-coordinates as intermediate
             * -> setting state manually in order to keep crystal valid
             */
            asFmt(AtomFmt::Crystal).getNat();
            cell->cellvec = vec;
            cell->invvec = inv;
            constexpr auto crystal = static_cast<size_t>(AtomFmt::Crystal);
            const auto target = static_cast<size_t>(at_fmt);
            atoms->coordinates[target] =
                    formatAll(atoms->coordinates[crystal],
                              AtomFmt::Crystal, at_fmt);
            for (size_t i=0; i<nAtFmt; ++i){
                if ((i == crystal) || (i == target)) continue;
                atoms->coord_outdated[i] = true;
            }
        }
    } else {
        if (at_fmt == AtomFmt::Crystal) {
            /*
             * Crystal needs to be reformatted
             *
             * Determine a valid buffer or use Bohr
             * -> setting state manually in order to keep valid
             */
            evaluateCache();
            size_t buf = nAtFmt;
            for (size_t i=0; i<nAtFmt; ++i){
                if (atoms->coord_outdated[i] == false){
                    buf = i;
                }
            }
            if (buf == nAtFmt) {
                buf = static_cast<size_t>(AtomFmt::Bohr);
                asFmt(AtomFmt::Bohr).getNat();
            }
            cell->cellvec = vec;
            cell->invvec = inv;
            constexpr auto target = static_cast<size_t>(AtomFmt::Crystal);
            atoms->coordinates[target] =
                    formatAll(atoms->coordinates[buf],
                              static_cast<AtomFmt>(buf), AtomFmt::Crystal);
            for (size_t i=0; i<nAtFmt; ++i){
                if ((i == target) || (i == buf)) continue;
                atoms->coord_outdated[i] = true;
            }
        } else {
            /*
             * All but crystal stay as-is
             * bonds should be reevaluated
             */
            cell->cellvec = vec;
            cell->invvec = inv;
            bonds->outdated = true;
        }
    }
}

void Step::setCellDim(float cdm, CdmFmt fmt, bool scale)
{
    if(!(cdm>0))throw Error("Step::setCellDim(): "
                            "cell-dimension must be positive");
    enableCell(true);
    float ratio = -1;
    if (scale && (at_fmt < AtomFmt::Crystal)) {
        ratio = cdm / getCellDim(fmt);
    } else if (!scale && (at_fmt >= AtomFmt::Crystal)) {
        ratio = getCellDim(fmt) / cdm;
    }
    if (ratio > 0) {
        evaluateCache();
        for(auto& c: atoms->coordinates[static_cast<size_t>(at_fmt)]){
            c *= ratio;
        }
        atoms->coord_changed[static_cast<size_t>(at_fmt)] = true;
    } else if (bonds->level == BondLevel::Cell) {
        bonds->outdated = true;
    }
    if (fmt == CdmFmt::Bohr) {
        cell->dimBohr = cdm;
        cell->dimAngstrom = cdm*bohrrad;
    } else {
        cell->dimAngstrom = cdm;
        cell->dimBohr = cdm*invbohr;
    }
}

Step& StepFormatter::asFmt(AtomFmt f)
{
    return step.asFmt(f);
}

const Step& StepFormatter::asFmt(AtomFmt f) const
{
    return step.asFmt(f);
}

Step& StepProper::asFmt(AtomFmt f)
{
    return this->*fmtmap[static_cast<size_t>(f)];
}

const Step& StepProper::asFmt(AtomFmt f) const
{
    return this->*fmtmap[static_cast<size_t>(f)];
}

void StepProper::setFmt(AtomFmt at_fmt, bool scale){
    if(at_fmt == this->at_fmt){ return; }
    auto source = static_cast<size_t>(this->at_fmt);
    auto target = static_cast<size_t>(at_fmt);
    AtomList& atoms = *this->atoms;
    if (scale) {
        /*
         * 'scaling' means we keep all coordinates as-is,
         * but change the fmt
         */
        evaluateCache(); // pull in changes from elsewhere
        std::swap(atoms.coordinates[source], atoms.coordinates[target]);
        atoms.coord_changed[target] = true; // make sure other caches are invalidated
    } else {
        /*
         * without scaling, we just change which buffer we're pulling data from
         * (after making sure it's up to date by
         *  calling getNat() which evaluates the cache)
         */
        asFmt(at_fmt).getNat();
    }
    if (at_fmt >= AtomFmt::Crystal) enableCell(true);
    this->at_fmt = at_fmt;
}

Step::Step(std::shared_ptr<PseMap> p, AtomFmt f, std::shared_ptr<BondList> b,
         std::shared_ptr<CellData> c, std::shared_ptr<std::string> s,
         std::shared_ptr<AtomList> a)
    : StepBase{p,f,b,c,s}, atoms{a} {}

Step::Step(const Step& s)
    : StepBase{s.pse, s.at_fmt,
               std::make_shared<BondList>(*s.bonds),
               std::make_shared<CellData>(*s.cell),
               std::make_shared<std::string>(*s.comment)},
      atoms{std::make_shared<AtomList>(*s.atoms)} {}

Step& Step::operator=(const Step& s)
{
    pse = s.pse;
    at_fmt = s.at_fmt;
    atoms = std::make_shared<AtomList>(*s.atoms);
    bonds = std::make_shared<BondList>(*s.bonds);
    cell = std::make_shared<CellData>(*s.cell);
    comment = std::make_shared<std::string>(*s.comment);
    return *this;
}

StepFormatter::StepFormatter(StepProper& step, AtomFmt at_fmt)
    : Step{step.pse, at_fmt, step.bonds, step.cell,
           step.comment, step.atoms}, step{step} {}

StepProper::StepProper(std::shared_ptr<PseMap> pse, AtomFmt at_fmt,
                       std::string comment)
    : Step{pse, at_fmt,
           std::make_shared<BondList>(),
           std::make_shared<CellData>(),
           std::make_shared<std::string>(comment),
           std::make_shared<AtomList>()} {}

StepProper::StepProper(const StepProper& s)
    : Step{s} {}

StepProper& StepProper::operator=(const StepProper& s)
{
    static_cast<Step&>(*this) = s;
    return *this;
}

constexpr StepFormatter StepProper::* StepProper::fmtmap[];

void StepFormatter::evaluateCache() const
{
    // let Step sort out overall situation
    step.evaluateCache();
    // if still outdated, pull in changes from Step
    const auto fmt = static_cast<size_t>(at_fmt);
    if (atoms->coord_outdated[fmt]) {
        atoms->coordinates[fmt] =
                formatAll(atoms->coordinates[static_cast<size_t>(step.at_fmt)],
                          step.at_fmt, at_fmt);
        atoms->coord_outdated[fmt] = false;
    }
}

void StepProper::evaluateCache() const
{
    auto fmt = static_cast<size_t>(this->at_fmt);
    auto& atoms = *this->atoms;
    // if there are (only) local changes, invalidate format-caches
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
        // if there are changes in (only one) formatter,
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
                                               static_cast<AtomFmt>(fmt));
            for (size_t i=0; i<nAtFmt; ++i){
                if((i == fmt) || (i == fmt_changed)){
                    atoms.coord_changed[i] = false;
                    atoms.coord_outdated[i] = false;
                }else{
                    atoms.coord_outdated[i] = true;
                }
            }
        }
        bonds->outdated = true;
    }
    // if some atom-types changed, update pse pointers
    if(atoms.prop_changed){
        if(size_t nat = atoms.names.size()){
            for(size_t i=0; i<nat; ++i){
                atoms.pse[i] = &(*pse)[atoms.names[i]];
            }
        }
        atoms.prop_changed = false;
        bonds->outdated = true;
    }
}
