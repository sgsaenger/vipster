#include "stepsel.h"
using namespace Vipster;

StepSelection::StepSelection(std::shared_ptr<PseMap> pse, AtomFmt fmt,
                             std::shared_ptr<BondList> bonds,
                             std::shared_ptr<CellData> cell,
                             std::shared_ptr<std::string> comment,
                             std::shared_ptr<AtomSelection> selection,
                             std::shared_ptr<Criteria> criteria,
                             Step& step)
    : StepBase{std::move(pse), fmt, std::move(bonds), std::move(cell), std::move(comment)},
      selection{selection},
      criteria{criteria},
      step{step}
{}

StepSelection::StepSelection(const StepSelection& s)
    : StepBase{s.pse, s.at_fmt,
               std::make_shared<BondList>(*s.bonds),
               std::make_shared<CellData>(*s.cell),
               std::make_shared<std::string>(*s.comment)},
      selection{std::make_shared<AtomSelection>(*s.selection)},
      criteria{std::make_shared<Criteria>(*s.criteria)},
      step{s.step}
{}

void StepSelection::evaluateCache() const
{
    step.evaluateCache();
    // TODO: evaluate selection
}

Criteria& StepSelection::getCriteria() noexcept
{
    return *criteria;
}

const Criteria& StepSelection::getCriteria() const noexcept
{
    return *criteria;
}

// Atoms
size_t StepSelection::getNat() const noexcept
{
    return selection->indices.size();
}

Atom StepSelection::operator[](size_t i)
{
    evaluateCache();
    return *iterator{selection, at_fmt, i};
}

StepSelection::iterator StepSelection::begin() noexcept
{
    evaluateCache();
    return iterator{selection, at_fmt, 0};
}

StepSelection::constIterator StepSelection::begin() const noexcept
{
    evaluateCache();
    return constIterator{selection, at_fmt, 0};
}

StepSelection::constIterator StepSelection::cbegin() const noexcept
{
    evaluateCache();
    return constIterator{selection, at_fmt, 0};
}

StepSelection::iterator StepSelection::end() noexcept
{
    evaluateCache();
    return iterator{selection, at_fmt, 0};
}

StepSelection::constIterator StepSelection::end() const noexcept
{
    evaluateCache();
    return constIterator{selection, at_fmt, 0};
}

StepSelection::constIterator StepSelection::cend() const noexcept
{
    evaluateCache();
    return constIterator{selection, at_fmt, 0};
}

/*
 * Formatter specific methods
 */
StepSelFormatter::StepSelFormatter(StepSelProper& selector, AtomFmt fmt)
    : StepSelection{selector.pse, fmt, selector.bonds, selector.cell,
                    selector.comment, selector.selection, selector.criteria,
                    selector.step.asFmt(fmt)},
      selector{selector}
{}

StepSelection& StepSelFormatter::asFmt(AtomFmt f)
{
    return selector.asFmt(f);
}

const StepSelection& StepSelFormatter::asFmt(AtomFmt f) const
{
    return selector.asFmt(f);
}

/*
 * Proper specific methods
 */
StepSelProper::StepSelProper(StepProper& step)
    : StepSelection{step.pse, step.getFmt(),
                    std::make_shared<BondList>(),
                    step.cell,
                    std::make_shared<std::string>(step.getComment()),
                    std::make_shared<AtomSelection>(),
                    std::make_shared<Criteria>(),
                    step}
{}

StepSelProper::StepSelProper(const StepSelProper& s)
    : StepSelection{s}
{}

StepSelection& StepSelProper::asFmt(AtomFmt f)
{
    return this->*fmtmap[static_cast<size_t>(f)];
}

const StepSelection& StepSelProper::asFmt(AtomFmt f) const
{
    return this->*fmtmap[static_cast<size_t>(f)];
}

constexpr StepSelFormatter StepSelProper::* StepSelProper::fmtmap[];
