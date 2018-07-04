#include "plugin.h"

using namespace Vipster;

IO::XYZConfig::XYZConfig(std::string n, Mode m, Data d)
    : BaseConfig{n}, filemode{m}, atomdata{d}
{}

std::unique_ptr<BaseConfig> IO::XYZConfig::copy()
{
    return std::make_unique<IO::XYZConfig>(*this);
}
