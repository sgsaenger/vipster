#include "vipsterapplication.h"

using namespace Vipster;

Application::Application()
{

}

void Application::setActiveStep(Step &step, Step::selection &sel)
{
    curStep = &step;
    curSel = &sel;
    emit activeStepChanged(step, sel);
}

void Application::newMol(Molecule &&mol){
    molecules.emplace_back(std::move(mol));
    molecules.back().getPTE().root = &config.periodicTable;

    emit molListChanged(molecules);
}

void Application::newPreset(Preset &&p){
    presets.emplace_back(std::move(p));

    emit presetListChanged();
}

void Application::newParameter(Parameter &&p){
    parameters.emplace_back(std::move(p));

    emit parameterListChanged();
}

void Application::newData(std::unique_ptr<const BaseData> &&d){
    data.emplace_back(std::move(d));

    emit dataListChanged();
}

void Application::newIOData(IOTuple &&t){
    newMol(std::move(std::get<Molecule>(t)));

    if (auto &param = std::get<std::optional<Parameter>>(t)) {
        newParameter(std::move(*param));
    }

    for (auto &extra: std::get<DataList>(t)) {
        newData(std::move(extra));
    }
}

void Application::selectionToCopy(){
    copyBuf = std::make_unique<Step::selection>(*curSel);
}
