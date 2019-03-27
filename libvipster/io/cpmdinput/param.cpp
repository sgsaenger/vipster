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

IOFmt IO::CPParam::getFmt() const
{
    return IOFmt::CPI;
}

std::unique_ptr<IO::BaseParam> IO::CPParam::copy() const
{
    return std::make_unique<IO::CPParam>(*this);
}

const std::vector<std::pair<std::string, IO::CPParam::Section IO::CPParam::*>> IO::CPParam::str2section{
    {"&INFO", &IO::CPParam::info},
    {"&CPMD", &IO::CPParam::cpmd},
    {"&SYSTEM", &IO::CPParam::system},
    {"&PIMD", &IO::CPParam::pimd},
    {"&PATH", &IO::CPParam::path},
    {"&PTDDFT", &IO::CPParam::ptddft},
    {"&ATOMS", &IO::CPParam::atoms},
    {"&DFT", &IO::CPParam::dft},
    {"&PROP", &IO::CPParam::prop},
    {"&RESP", &IO::CPParam::resp},
    {"&LINRES", &IO::CPParam::linres},
    {"&TDDFT", &IO::CPParam::tddft},
    {"&HARDNESS", &IO::CPParam::hardness},
    {"&CLASSIC", &IO::CPParam::classic},
    {"&EXTE", &IO::CPParam::exte},
    {"&VDW", &IO::CPParam::vdw},
    {"&QMMM", &IO::CPParam::qmmm},
};

void IO::CPParam::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::CPParam::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}

void IO::from_json(const nlohmann::json& j, IO::CPParam& p)
{
    for(const auto& pair: CPParam::str2section){
        p.*pair.second = j.value(pair.first, CPParam::Section{});
    }
}

void IO::to_json(nlohmann::json& j,const CPParam& p)
{
    for(const auto& pair: CPParam::str2section){
        j[pair.first] = p.*pair.second;
    }
}
