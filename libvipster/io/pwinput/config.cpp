#include "config.h"

using namespace Vipster;

IO::PWConfig::PWConfig(std::string name, PWConfig::AtomFmt atoms, CellFmt cell)
    : BaseConfig{name}, atoms{atoms}, cell{cell}
{}

std::unique_ptr<BaseConfig> IO::PWConfig::copy()
{
    return std::make_unique<IO::PWConfig>(*this);
}
