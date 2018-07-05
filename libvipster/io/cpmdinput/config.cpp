#include "config.h"

using namespace Vipster;

IO::CPConfig::CPConfig(std::string name, bool angstrom, Scale scale)
    : BaseConfig{name}, angstrom{angstrom}, scale{scale}
{}

std::unique_ptr<BaseConfig> IO::CPConfig::copy()
{
    return std::make_unique<IO::CPConfig>(*this);
}
