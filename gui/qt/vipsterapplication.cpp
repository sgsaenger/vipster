#include "vipsterapplication.h"

using namespace Vipster;

Application::Application()
{
    // Ensure bonds and overlap information is up to date
    auto updateBonds = [&](Step &step){
        if (getState(step).automatic_bonds ||
            config_val.settings.overlap.val) {
            step.generateBonds(!getState(step).automatic_bonds);
        }
    };
    connect(this, &Application::stepChanged, updateBonds);
}

const ConfigState& Application::config()
{
    return config_val;
}

Molecule& Application::curMol()
{
    return *pCurMol;
}

void Application::setActiveMol(Molecule &m)
{
    pCurMol = &m;

    emit activeMolChanged(m);
}

const Step& Application::curStep()
{
    return *pCurStep;
}

const Step::selection& Application::curSel()
{
    return *pCurSel;
}

void Application::setActiveStep(Step &step, Step::selection &sel)
{
    pCurStep = &step;
    pCurSel = &sel;

    getState(step);
    emit activeStepChanged(step, sel);
}

Application::StepState& Application::getState(const Step &step)
{
    // initialize state if required
    if (stepdata.find(&step) == stepdata.end()) {
        stepdata.emplace(&step, StepState{
            // when bonds are present (e.g. in file), use manual mode
            step.getBonds().empty(),
            // default formatter is passthrough
            const_cast<Step&>(step).asFmt(step.getFmt())
        });
    }
    return stepdata.at(&step);
}

void Application::newMol(Molecule &&mol){
    molecules.emplace_back(std::move(mol));
    molecules.back().getPTE().root = &config_val.periodicTable;

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
    copyBuf = std::make_unique<Step::selection>(*pCurSel);
    emit copyBufChanged(*copyBuf);
}
