#include "param.h"
#include "plugin.h"

using namespace Vipster;

IO::CPParam::CPParam(Section info, Section cpmd, Section system, Section pimd,
                     Section path, Section ptddft, Section atoms, Section dft,
                     Section prop, Section resp, Section linres, Section tddft,
                     Section hardness, Section classic, Section exte, Section vdw,
                     Section qmmm, std::string PPPrefix, std::string PPSuffix,
                     std::string PPNonlocality)
    : info{info}, cpmd{cpmd}, system{system}, pimd{pimd}, path{path},
      ptddft{ptddft}, atoms{atoms}, dft{dft}, prop{prop}, resp{resp}, linres{linres},
      tddft{tddft}, hardness{hardness}, classic{classic}, exte{exte}, vdw{vdw}, qmmm{qmmm},
      PPPrefix{PPPrefix}, PPSuffix{PPSuffix}, PPNonlocality{PPNonlocality}
{}

const IO::Plugin *IO::CPParam::getFmt() const
{
    return &CPInput;
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

void IO::CPParam::parseJson(const nlohmann::json& j)
{
    for(const auto& pair: CPParam::str2section){
        this->*pair.second = j.value(pair.first, CPParam::Section{});
    }
    PPPrefix = j.value("PPPrefix", "");
    PPSuffix = j.value("PPSuffix", "");
    PPNonlocality = j.value("PPNonlocality", "");
}

nlohmann::json IO::CPParam::toJson() const
{
    nlohmann::json j;
    for(const auto& pair: CPParam::str2section){
        j[pair.first] = this->*pair.second;
    }
    j["PPPrefix"] = PPPrefix;
    j["PPSuffix"] = PPSuffix;
    j["PPNonlocality"] = PPNonlocality;
    return j;
}
