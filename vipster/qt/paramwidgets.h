#ifndef PARAMWIDGETS_H
#define PARAMWIDGETS_H

#include "paramwidgets/pwparam.h"
#include "paramwidgets/orcaparam.h"
#include "paramwidgets/cpparam.h"

inline std::map<Vipster::IOFmt, ParamBase*> makeParamWidgets()
{
    return {
        {Vipster::IOFmt::PWI, new PWParam()},
        {Vipster::IOFmt::ORCA, new ORCAParam()},
        {Vipster::IOFmt::CPI, new CPParam()},
    };
}

#endif // PARAMWIDGETS_H
