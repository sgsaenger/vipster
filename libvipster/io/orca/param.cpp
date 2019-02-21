#include "param.h"

using namespace Vipster;

IO::OrcaParam::OrcaParam(Header header)
    : header{header}
{}

IOFmt IO::OrcaParam::getFmt() const
{
    return IOFmt::ORCA;
}

std::unique_ptr<IO::BaseParam> IO::OrcaParam::copy() const
{
    return std::make_unique<IO::OrcaParam>(*this);
}

void IO::OrcaParam::parseJson(const nlohmann::json& j)
{
    j.get_to(header);
}

nlohmann::json IO::OrcaParam::toJson() const
{
    nlohmann::json j;
    j = header;
    return j;
}
