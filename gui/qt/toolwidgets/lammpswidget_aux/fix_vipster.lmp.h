#ifndef LMP_FIX_VIPSTER_H
#define LMP_FIX_VIPSTER_H

#include "fix.h"

#include "vipster/molecule.h"

namespace LAMMPS_NS {

Fix* mkFixVipster(class LAMMPS*, int narg, char **arg);

class FixVipster: public Fix {
public:
    FixVipster(class LAMMPS*, int narg, char **arg);
    int setmask() override;
    void post_run() override;
    void min_post_force(int) override;
    void final_integrate() override;
    MPI_Comm intercomm{MPI_COMM_NULL};
    Vipster::Molecule *mol{nullptr};
private:
    void sendStep();
    int nevery{1};
};

}

#endif // LMP_FIX_VIPSTER_H
