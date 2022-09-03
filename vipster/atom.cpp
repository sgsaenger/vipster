#include "atom.h"

using namespace Vipster;

bool Vipster::atomFmtRelative(AtomFmt f)
{
    return (f == AtomFmt::Crystal) || (f == AtomFmt::Alat);
}

bool Vipster::atomFmtAbsolute(AtomFmt f)
{
    return (f > AtomFmt::Alat) && (f < static_cast<AtomFmt>(detail::toAngstrom.size()));
}

bool Vipster::AtomProperties::operator==(const AtomProperties &p) const
{
    return std::tie(charge, flags, forces)
           ==
           std::tie(p.charge, p.flags, p.forces);
}

detail::CoordConverter Vipster::detail::makeConverter(const AtomContext &source,
                                                      const AtomContext &target)
{
    switch(source.fmt){
    case AtomFmt::Crystal:
        switch(target.fmt){
        case AtomFmt::Crystal:
            if(source.cell == target.cell){
                return [](const Vec &v){return v;};
            }else{
                return [&](const Vec &v){return v * (source.cell->matrix * source.cell->dimension)
                                                  * (target.cell->inverse / target.cell->dimension);};
            }
        case AtomFmt::Alat:
            return [&](const Vec &v){return v * (source.cell->matrix * source.cell->dimension)
                                              / target.cell->dimension;};
        default:
            return [&](const Vec &v){return v * (source.cell->matrix * source.cell->dimension)
                                              * detail::fromAngstrom[static_cast<size_t>(target.fmt)];};
        }
    case AtomFmt::Alat:
        switch(target.fmt){
        case AtomFmt::Crystal:
            if(source.cell->dimension == target.cell->dimension){
                return [&](const Vec &v){return v * target.cell->inverse;};
            }else{
                return [&](const Vec &v){return v * source.cell->dimension
                                                  * (target.cell->inverse / target.cell->dimension);};
            }
        case AtomFmt::Alat:
            if(source.cell->dimension == target.cell->dimension){
                return [](const Vec &v){return v;};
            }else{
                return [&](const Vec &v){return v * source.cell->dimension / target.cell->dimension;};
            }
        default:
            return [&](const Vec &v){return v * source.cell->dimension
                                              * detail::fromAngstrom[static_cast<size_t>(target.fmt)];};
        }
    default: // absolute coordinates
        switch(target.fmt){
        case AtomFmt::Crystal:
            return [&](const Vec &v){return v * detail::toAngstrom[static_cast<size_t>(source.fmt)]
                                              * (target.cell->inverse / target.cell->dimension);};
        case AtomFmt::Alat:
            return [&](const Vec &v){return v * detail::toAngstrom[static_cast<size_t>(source.fmt)]
                                              / target.cell->dimension;};
        default:
            if(source.fmt == target.fmt){
                return [](const Vec &v){return v;};
            }else{
                return [&](const Vec &v){return v * detail::toAngstrom[static_cast<size_t>(source.fmt)]
                                                  * detail::fromAngstrom[static_cast<size_t>(target.fmt)];};
            }
        }
    }
}
