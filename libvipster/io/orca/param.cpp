#include "param.h"

using namespace Vipster;

IO::OrcaParam::OrcaParam(std::string name, Header header)
    : BaseParam{name}, header{header}
{}

IOFmt IO::OrcaParam::getFmt() const
{
    return IOFmt::ORCA;
}

std::unique_ptr<IO::BaseParam> IO::OrcaParam::copy() const
{
    return std::make_unique<IO::OrcaParam>(*this);
}

void IO::to_json(nlohmann::json& j, const OrcaParam& p)
{
    j = p.header;
}

void IO::from_json(const nlohmann::json& j, OrcaParam& p)
{
    j.get_to(p.header);
}

void IO::OrcaParam::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::OrcaParam::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}
