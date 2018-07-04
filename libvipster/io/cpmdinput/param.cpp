#include "param.h"

using namespace Vipster;

IO::CPParam::CPParam(std::string name, Section info, Section cpmd, Section system,
                     Section pimd, Section path, Section ptddft, Section atoms,
                     Section dft, Section prop, Section resp, Section linres,
                     Section tddft, Section hardness, Section classic, Section vdw, Section qmmm)
    : BaseParam{name}, info{info}, cpmd{cpmd}, system{system}, pimd{pimd}, path{path},
      ptddft{ptddft}, atoms{atoms}, dft{dft}, prop{prop}, resp{resp}, linres{linres},
      tddft{tddft}, hardness{hardness}, classic{classic}, vdw{vdw}, qmmm{qmmm}
{}

std::unique_ptr<BaseParam> IO::CPParam::copy()
{
    return std::make_unique<IO::CPParam>(*this);
}
