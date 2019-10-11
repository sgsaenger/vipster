#ifndef MAINWIDGETS_H
#define MAINWIDGETS_H

#include "mainwidgets/molwidget.h"
#include "mainwidgets/paramwidget.h"
#include "mainwidgets/configwidget.h"
#include "mainwidgets/periodictablewidget.h"
#include "mainwidgets/settingswidget.h"
#include "mainwidgets/datawidget.h"

inline std::vector<std::pair<BaseWidget*, QString>> makeMainWidgets(QWidget* parent)
{
    return {
        {new MolWidget(parent), "Molecule"},
        {new ParamWidget(parent), "Parameter"},
        {new ConfigWidget(parent), "Config"},
        {new SettingsWidget(parent), "Settings"},
        {new PeriodicTableWidget(parent), "Periodic Table (Molecule)"},
        {new PeriodicTableWidget(parent, true), "Periodic Table (Global)"},
        {new DataWidget(parent), "Additional Data"},
    };
}

#endif // MAINWIDGETS_H
