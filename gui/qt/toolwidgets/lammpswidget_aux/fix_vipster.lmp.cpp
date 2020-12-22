#include "viewport.h"
#include "fix_vipster.lmp.h"

#include "lammps.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "update.h"
#include "utils.h"

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
    nevery = utils::inumeric(FLERR, arg[3], false, lmp);
    if (nevery < 0) error->all(FLERR, "Illegal fix vipster command");
}

int FixVipster::setmask()
{
    return FINAL_INTEGRATE | MIN_POST_FORCE | POST_RUN;
}

void FixVipster::post_run()
{
    // TODO
    // MIN: report convergence?
    // MD: send last step here instead of final_integrate?
}

void FixVipster::min_post_force(int)
{
    sendStep();
}

void FixVipster::final_integrate()
{
    if(!(update->ntimestep % nevery) || update->ntimestep == update->laststep){
        sendStep();
    }
}

void FixVipster::sendStep()
{
    if(comm->nprocs > 1 && intercomm == MPI_COMM_NULL){
        throw LAMMPSException{"Fix vipster run in parallel and not setup for MPI reporting."};
    }else if(comm->nprocs == 1 && mol == nullptr){
        throw LAMMPSException{"Fix vipster run in serial and not setup for in-memory reporting."};
    }
    if(mol){
        auto &step = mol->newStep(mol->getStep(mol->getNstep()-1));
        if(atom->natoms != step.getNat()){
            throw LAMMPSException{"Number of atoms differs between LAMMPS and Vipster. Aborting."};
        }
        for(size_t i=0; i<atom->natoms; ++i){
            auto at = step.at(atom->tag[i]-1);
            at.coord[0] = atom->x[i][0];
            at.coord[1] = atom->x[i][1];
            at.coord[2] = atom->x[i][2];
            at.properties->forces[0] = atom->x[i][0];
            at.properties->forces[1] = atom->x[i][1];
            at.properties->forces[2] = atom->x[i][2];
        }
    }
    if(intercomm != MPI_COMM_NULL){
        // TODO
        // me==0 send status over intercomm -> newStep or done
        // how is "done" defined, especially in MIN?
        // MPI_Barrier on WORLD?
        // then MPI_Gather x and f over intercomm
    }
}
