#include "forcefield.lmp.h"
#include "uff.lmp.h"

ForceFields defaultForcefields()
{
    return {{"UFF", &UFF}};
};
