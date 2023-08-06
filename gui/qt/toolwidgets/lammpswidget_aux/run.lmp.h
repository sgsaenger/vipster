#ifndef RUN_LMP_H
#define RUN_LMP_H

#include "vipster/molecule.h"

#include <vector>
#include <string>

#include <mpi.h> // can be LAMMPS' dummy implementation

namespace Vipster::Lammps{
struct runParams{
    enum class Mode{Min, MD};
    Mode mode;
    // Parallelization
    int MPI{1}, OMP{0}, GPU{0};
    // MD + Min
    size_t nstep;
    // Min
    size_t neval{0};
    double etol{0}, ftol{0};
};

std::pair<int, std::string> runMaster(const std::string &dir, runParams params, Molecule *mol);
void runSlave();
void run(const std::string &dir, runParams params, MPI_Comm intercomm, Molecule *mol=nullptr);
}

#endif // RUN_LMP_H
