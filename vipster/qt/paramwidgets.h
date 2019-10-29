#ifndef PARAMWIDGETS_H
#define PARAMWIDGETS_H

#include "paramwidgets/pwparam.h"
#include "paramwidgets/orcaparam.h"
#include "paramwidgets/cpparam.h"

inline std::map<const Vipster::IO::Plugin*, ParamBase*> makeParamWidgets()
{
    return {
        {&Vipster::IO::PWInput, new PWParam()},
        {&Vipster::IO::OrcaInput, new ORCAParam()},
        {&Vipster::IO::CPInput, new CPParam()},
    };
}

#endif // PARAMWIDGETS_H
