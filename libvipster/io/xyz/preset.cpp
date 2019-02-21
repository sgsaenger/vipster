#include "preset.h"

using namespace Vipster;

IO::XYZPreset::XYZPreset(Mode m, Data d)
    : filemode{m}, atomdata{d}
{}

IOFmt IO::XYZPreset::getFmt() const
{
    return IOFmt::XYZ;
}

std::unique_ptr<IO::BasePreset> IO::XYZPreset::copy() const
{
    return std::make_unique<IO::XYZPreset>(*this);
}

void IO::XYZPreset::parseJson(const nlohmann::json& j)
{
    filemode = j.value("filemode", IO::XYZPreset::Mode::Step);
    atomdata = j.value("atomdata", IO::XYZPreset::Data::None);
}

nlohmann::json IO::XYZPreset::toJson() const
{
    nlohmann::json j;
    j["filemode"] = filemode;
    j["atomdata"] = atomdata;
    return j;
}
