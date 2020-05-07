#ifndef TOOLS_H
#define TOOLS_H

#include "basewidget.h"

#include "toolwidgets/cellmodwidget.h"
#include "toolwidgets/definewidget.h"
#include "toolwidgets/millerwidget.h"
#include "toolwidgets/pickwidget.h"
#include "toolwidgets/pinwidget.h"
#include "toolwidgets/scriptwidget.h"
#ifdef USE_PYTHON
#include "toolwidgets/pythonwidget.py.h"
#endif
#ifdef USE_LAMMPS
#include "toolwidgets/lammpswidget.lmp.h"
#endif

inline std::vector<std::pair<BaseWidget*, QString>> makeToolWidgets(QWidget* parent)
{
    return {
        {new PickWidget(parent), "Selected Atoms"},
        {new DefineWidget(parent), "Filter Atoms"},
        {new MillerWidget(parent), "Lattice Planes"},
        {new PinWidget(parent), "Pin Steps"},
        {new CellModWidget(parent), "Modify Cell"},
        {new ScriptWidget(parent), "Script"},
#ifdef USE_PYTHON
        {new PythonWidget(parent), "Python"},
#endif
#ifdef USE_LAMMPS
        {new LammpsWidget(parent), "LAMMPS"},
#endif
    };
}

#endif // TOOLS_H
