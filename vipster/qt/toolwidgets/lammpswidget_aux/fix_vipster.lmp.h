#ifndef LMP_FIX_VIPSTER_H
#define LMP_FIX_VIPSTER_H

#include "lammps/fix.h"
#include "mainwindow.h"

namespace LAMMPS_NS {

Fix* mkFixVipster(class LAMMPS*, int narg, char **arg);

class FixVipster: public Fix {
public:
    FixVipster(class LAMMPS*, int narg, char **arg);
    ~FixVipster();
    int setmask() override;
    void min_post_force(int) override;
    void final_integrate() override;
    void init_vipster(MainWindow *mw, const std::string &name);
private:
    MainWindow *master{nullptr};
    Vipster::Molecule *molecule{nullptr};
    int freq{1};
    int mol_idx{-1};
};

}

#endif // LMP_FIX_VIPSTER_H
