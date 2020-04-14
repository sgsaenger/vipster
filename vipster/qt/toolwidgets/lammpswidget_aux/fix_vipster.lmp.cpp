#include "viewport.h"
#include "fix_vipster.lmp.h"

#include "lammps/lammps.h"
#include "lammps/atom.h"
#include "lammps/error.h"
#include "lammps/force.h"
#include "lammps/update.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace Vipster;

Fix* LAMMPS_NS::mkFixVipster(LAMMPS *lmp, int narg, char **arg)
{
    return new FixVipster{lmp, narg, arg};
}

FixVipster::FixVipster(LAMMPS *lmp, int narg, char **arg):
    Fix{lmp, narg, arg}
{
    if (narg != 4) error->all(FLERR, "Illegal fix vipster command");
    nevery = force->inumeric(FLERR, arg[3]);
    if (nevery < 0) error->all(FLERR, "Illegal fix vipster command");
}
FixVipster::~FixVipster()
{
    master->curVP->moldata[molecule].curStep = molecule->getNstep();
    master->curVP->setMol(master->molecules.size()-1);
}

int FixVipster::setmask()
{
    return FINAL_INTEGRATE | MIN_POST_FORCE;
}

void FixVipster::min_post_force(int)
{
    final_integrate();
}

void FixVipster::final_integrate()
{
    // TODO: update view when executed in thread
    if (!master || !molecule)
        throw LAMMPSException{"fix vipster not correctly initialized by Vipster."};
    // transfer data on requested steps
    if(!(update->ntimestep % nevery) || update->ntimestep == update->laststep){
        copyCurStep();
    }
}

void FixVipster::copyCurStep()
{
    // copy latest step
    auto &step = molecule->newStep(molecule->getStep(molecule->getNstep()-1));
    // collect x and f
    if (atom->natoms != step.getNat())
        throw LAMMPSException{"Number of atoms differs between LAMMPS and Vipster. Aborting."};
    for(bigint i=0; i < atom->natoms; ++i){
        auto at = step.at(atom->tag[i]-1);
        at.coord[0] = atom->x[i][0];
        at.coord[1] = atom->x[i][1];
        at.coord[2] = atom->x[i][2];
        at.properties->forces[0] = atom->x[i][0];
        at.properties->forces[1] = atom->x[i][1];
        at.properties->forces[2] = atom->x[i][2];
    }
}

void FixVipster::init_vipster(MainWindow *mw, const std::string &name)
{
    master = mw;
    // create new trajectory
    master->newMol(Vipster::Molecule{*master->curStep, master->curMol->name + ' ' + name});
    // save pointer of own trajectory
    molecule = master->curMol;
}
