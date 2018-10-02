#include "step.h"
#include "stepsel.h"

using namespace Vipster;

void Step::enableCell(bool b) noexcept
{
    //TODO: create minimum containing cell?
    cell->enabled = b;
}

size_t Step::getNat() const noexcept{
    return atoms->names.size();
}

void Step::newAtom(){
    AtomList& al = *atoms;
    // Coordinates
    al.coordinates[static_cast<size_t>(at_fmt)].emplace_back();
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.names.emplace_back();
    al.name_changed = true;
    al.pse.push_back(&(*pse)[""]);
    // Properties
    al.properties.emplace_back();
    al.prop_changed = true;
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

void Step::newAtoms(const AtomList& atoms){
    size_t nat = getNat() + atoms.names.size();
    // Coordinates
    AtomList& al = *this->atoms;
    al.coordinates[static_cast<size_t>(at_fmt)].reserve(nat);
    al.coordinates[static_cast<size_t>(at_fmt)].insert(
                al.coordinates[static_cast<size_t>(at_fmt)].end(),
                atoms.coordinates[static_cast<size_t>(at_fmt)].begin(),
                atoms.coordinates[static_cast<size_t>(at_fmt)].end());
    al.coord_changed[static_cast<size_t>(at_fmt)] = true;
    // Type
    al.names.reserve(nat);
    al.names.insert(al.names.end(), atoms.names.begin(), atoms.names.end());
    al.name_changed = true;
    al.pse.resize(nat);
    // Properties
    al.properties.reserve(nat);
    al.properties.insert(al.properties.end(), atoms.properties.begin(), atoms.properties.end());
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

Step::iterator Step::begin() noexcept {
    return {atoms, at_fmt, 0};
}

Step::constIterator Step::begin() const noexcept {
    return {atoms, at_fmt, 0};
}

Step::constIterator Step::cbegin() const noexcept {
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

const AtomList& Step::getAtoms() const noexcept {
    return *atoms;
}

const std::vector<Vec>& Step::getCoords() const noexcept {
    return atoms->coordinates[static_cast<size_t>(at_fmt)];
}

const std::vector<std::string>& Step::getNames() const noexcept {
    return atoms->names;
}

const std::vector<PseEntry*>& Step::getPseEntries() const noexcept {
    return atoms->pse;
}

const std::vector<AtomProperties>& Step::getProperties() const noexcept {
    return atoms->properties;
}

Atom Step::operator[](size_t i) {
    return *iterator{atoms, at_fmt, i};
}

StepSelection Step::select(std::string filter)
{
    return StepSelection{*this, filter};
}

StepSelection Step::select(SelectionFilter filter)
{
    return StepSelection{*this, filter};
}

StepSelConst Step::select(std::string filter) const
{
    return StepSelConst{*this, filter};
}

StepSelConst Step::select(SelectionFilter filter) const
{
    return StepSelConst{*this, filter};
}

void Step::setCellVec(const Mat &vec, bool scale)
{
    Mat inv = Mat_inv(vec);
    cell->enabled = true;
    if (scale) {
        if (at_fmt == AtomFmt::Crystal) {
            /*
             * Crystal coordinates stay as-is
             */
            cell->cellvec = vec;
            cell->invvec = inv;
            atoms->coord_changed[static_cast<size_t>(at_fmt)] = true;
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
            for(size_t i=0; i<nAtFmt; ++i){
                if ((i == target) || (i == crystal)){
                    atoms->coord_changed[i] = false;
                    atoms->coord_outdated[i] = false;
                }else{
                    atoms->coord_outdated[i] = true;
                }
            }
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
            for(size_t i=0; i<nAtFmt; ++i){
                if ((i == target) || (i == buf)){
                    atoms->coord_changed[i] = false;
                    atoms->coord_outdated[i] = false;
                }else{
                    atoms->coord_outdated[i] = true;
                }
            }
        } else {
            /*
             * All but crystal stay as-is
             */
            cell->cellvec = vec;
            cell->invvec = inv;
            atoms->coord_changed[static_cast<size_t>(at_fmt)] = true;
        }
    }
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
        // make sure cache state is coherent
        atoms.coord_outdated[target] = false;
        atoms.coord_changed[target] = true;
    } else {
        /*
         * without scaling, we just change which buffer we're pulling data from
         */
        asFmt(at_fmt).evaluateCache();
    }
    if (at_fmt >= AtomFmt::Crystal){
        cell->enabled = true;
    }
    this->at_fmt = at_fmt;
}

Step::Step(std::shared_ptr<PseMap> p, AtomFmt f, std::shared_ptr<BondList> b,
         std::shared_ptr<CellData> c, std::shared_ptr<std::string> s,
         std::shared_ptr<AtomList> a)
    : StepBase{std::move(p),f,std::move(b),std::move(c),std::move(s)},
      atoms{std::move(a)} {}

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
    *atoms = *s.atoms;
    *bonds = *s.bonds;
    *cell = *s.cell;
    *comment = *s.comment;
    return *this;
}

StepFormatter::StepFormatter(StepProper& step, AtomFmt at_fmt)
    : Step{step.pse, at_fmt, step.bonds, step.cell,
           step.comment, step.atoms}, step{step} {}

StepProper::StepProper(std::shared_ptr<PseMap> pse, AtomFmt at_fmt,
                       std::string comment)
    : Step{std::move(pse), at_fmt,
           std::make_shared<BondList>(),
           std::make_shared<CellData>(),
           std::make_shared<std::string>(comment),
           std::make_shared<AtomList>()} {}

StepProper::StepProper(const StepProper& s)
    : Step{s} {}

StepProper& StepProper::operator=(const StepProper& s)
{
    *static_cast<Step*>(this) = s;
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
            bonds->outdated = true;
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
