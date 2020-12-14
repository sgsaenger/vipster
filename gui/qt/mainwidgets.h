#ifndef MAINWIDGETS_H
#define MAINWIDGETS_H

#include "mainwidgets/molwidget.h"
#include "mainwidgets/paramwidget.h"
#include "mainwidgets/presetwidget.h"
#include "mainwidgets/periodictablewidget.h"
#include "mainwidgets/settingswidget.h"
#include "mainwidgets/datawidget.h"

inline std::vector<std::pair<BaseWidget*, QString>> makeMainWidgets(QWidget* parent)
{
    return {
        {new MolWidget(parent), "Molecule"},
        {new ParamWidget(parent), "Parameter"},
        {new PresetWidget(parent), "IO-Presets"},
        {new SettingsWidget(parent), "Settings"},
        {new PeriodicTableWidget(parent, true), "Periodic Table"},
        {new DataWidget(parent), "Additional Data"},
    };
}

#endif // MAINWIDGETS_H
