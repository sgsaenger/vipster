#include "lammpswidget.lmp.h"
#include "ui_lammpswidget.lmp.h"
#include "../mainwindow.h"

#include <QMessageBox>

#ifndef LAMMPS_EXCEPTIONS
#define LAMMPS_EXCEPTIONS
#endif
#include "lammps.h"
#include "domain.h"
#include "error.h"
#include "exceptions.h"
#include "input.h"
#include "update.h"

using namespace Vipster;
using namespace LAMMPS_NS;

LammpsWidget::LammpsWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::LammpsWidget)
{
    ui->setupUi(this);
    // set validators for minimize settings
    ui->etolInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->ftolInput->setValidator(new QDoubleValidator{0, 100, 5});
    ui->maxItInput->setValidator(new QIntValidator{1, 10000000});
    ui->maxevalInput->setValidator(new QIntValidator{1, 10000000});
    // enable parallelization as available
#ifdef USE_MPI
    ui->MPISpin->setDisabled(true);
#else
    ui->MPISpin->setEnabled(true);
#endif
    ui->OMPSpin->setDisabled(true);
    ui->GPUSpin->setDisabled(true);
    int i=0;
    while(LAMMPS::installed_packages[i] != nullptr){
        if(std::string{"USER-OMP"} == LAMMPS::installed_packages[i]){
            ui->OMPSpin->setEnabled(true);
        }else if(std::string{"GPU"} == LAMMPS::installed_packages[i]){
            ui->GPUSpin->setEnabled(true);
        }
        ++i;
    }
//#define PAIR_CLASS
//#define PairStyle(key,Class) BARGL
//#include "style_pair.h"
//#undef PAIR_CLASS
}

LammpsWidget::~LammpsWidget()
{
    delete ui;
}

void LammpsWidget::on_runButton_clicked()
{
    std::array<char*, 5> lmparg{
        nullptr,
        "-screen", "none",
        "-log", "none"
    };
    LAMMPS lmp{5, lmparg.data(), MPI_COMM_WORLD};
    try{
        auto &curStep = *master->curStep;
        // initialize lammps
        lmp.init();
        // always use real units -> angstrom
        lmp.update->set_units("real");
        // create box
        const auto& vec = curStep.getCellVec();
        lmp.domain->boxlo[0] = 0;
        lmp.domain->boxlo[1] = 0;
        lmp.domain->boxlo[2] = 0;
        lmp.domain->boxhi[0] = vec[0][0];
        lmp.domain->boxhi[1] = vec[1][1];
        lmp.domain->boxhi[2] = vec[2][2];
        lmp.domain->xy = vec[1][0];
        lmp.domain->xz = vec[2][0];
        lmp.domain->yz = vec[2][1];
        lmp.domain->reset_box();
    }catch(LAMMPSAbortException &e){
        QMessageBox::critical(this, "Error in LAMMPS-Run", QString::fromStdString(e.message));
    }catch(LAMMPSException &e){
        QMessageBox::warning(this, "Error in LAMMPS-Run", QString::fromStdString(e.message));
    }
}

void LammpsWidget::on_helpButton_clicked()
{
    // TODO: explain this widget
}

std::vector<std::tuple<size_t, size_t, int>> aromaticity_criteria{
    // atomic number, coordination, delta ring-size (starting at 6)
    // boron group:
    {5,  3, +1},
    // carbon group:
    {6,  3,  0},
    {14, 3,  0},
    // nitrogen group:
    {7,  2,  0},
    {7,  3, -1},
    {15, 2,  0},
    {15, 3, -1},
    // oxygen group:
    {8,  2, -1},
    {16, 2, -1},
};

void LammpsWidget::on_ffPrepare_clicked()
{
    // TODO: UFF only atm, create abstraction
    const auto& origin = *master->curStep;
    Step local = origin;
    /* TODO:
     * 1. create adjacency list (std::map<std::vector<pointer to map-entries>>)
     * 2. search for aromatic rings (one combined step with dynamically adjusted ring-size?)
     * 3. assign atom-types based on coordination and aromaticity (yes/no/ignore)
     */
    using adjacency_node = std::vector<void*>;
    auto adjacency_list = [](const Step &step)
    {
        std::vector<adjacency_node> coord(step.getNat());
        for(const auto& b: step.getBonds()){
            coord[b.at1].push_back(&coord[b.at2]);
            coord[b.at2].push_back(&coord[b.at1]);
        }
        return coord;
    }(origin);
    auto search_aromatic_rings = [](auto& al, const Step& step){
        size_t target_size = 6;
        size_t current_size = 0;
        adjacency_node *const baseNode = &*al.begin();
        std::vector<adjacency_node*> indirect_al(al.size());
        for(size_t i=0; i<al.size(); ++i){
            indirect_al[i] = &al[i];
        }
        adjacency_node **curNode = &*indirect_al.begin();
        adjacency_node **curEnd = curNode + indirect_al.size();
        std::vector<std::tuple<adjacency_node**, adjacency_node**, int>> nodeStack;
        std::set<size_t> found_atoms;
        // iterate over all nodes
        while(true){
            const size_t curIdx = *curNode - baseNode;
            const std::string& curName = step[curIdx].name;
            // check if the current atom could possibly be aromatic
            const auto pos = std::find_if(aromaticity_criteria.begin(), aromaticity_criteria.end(), [&](const auto& type){
                return (std::get<0>(type) == step[*curNode - baseNode].type->Z)
                    && (std::get<1>(type) == (*curNode)->size());
            });
            /* abort-criteria:
             * - non aromatic atom
             * - aromatic atom, but size > target
             * - aromatic atom, size matches, but current atom != first atom
             * - aromatic atom, size too small, current atom already encountered
             *
             * success: aromatic atom, size matches, current atom == first atom
             * continue: aromatic atom, size too small, atom not encountered previously
             */
            if(pos != aromaticity_criteria.end()){
                if((current_size == target_size) &&
                   (curNode == std::get<0>(nodeStack[0]))){
                    // found an aromatic ring, success
                    for(const auto &[node, t1, t2]: nodeStack){
                        found_atoms.insert(*node - baseNode);
                    }
                }else if(current_size < target_size){
                    // too small
                    const auto tmp = std::find_if(nodeStack.begin(), nodeStack.end(), [&](const auto& tuple){return std::get<0>(tuple) == curNode;}); // FIXME: does not work as intended
                    if(tmp == nodeStack.end()){
                        // is not ring, increase depth
                        current_size += 1;
                        target_size += std::get<2>(*pos);
                        nodeStack.emplace_back(curNode, curEnd, std::get<2>(*pos));
                        auto tmpSize = (*curNode)->size();
                        curNode = reinterpret_cast<adjacency_node**>(&*(*curNode)->begin());
                        curEnd = curNode + tmpSize;
                        continue;
                    }
                }
            }
            // query next atom on current or previous depth
            while(++curNode == curEnd){
                // decrease depth
                current_size -= 1;
                int tmp;
                if(nodeStack.empty()){
                    // exhausted adjacency list, return results
                    return found_atoms;
                }
                std::tie(curNode, curEnd, tmp) = nodeStack.back();
                nodeStack.pop_back();
                target_size -= tmp;
            }
        }
    };//(adjacency_list);
    search_aromatic_rings(adjacency_list, origin);
}

std::vector<std::tuple<std::string, double, double, double, double, double, double>> UFF_Parameters{
{"H_",    0.354,  180.0, 2.886, 0.044,   12.0, 0.712},
{"H_b",   0.460,   83.5, 2.886, 0.044,   12.0, 0.712},
{"He4+4", 0.849,   90.0, 2.362, 0.056,  15.24, 0.098},
{"Li",    1.336,  180.0, 2.451, 0.025,   12.0, 1.026},
{"Be3+2", 1.074, 109.47, 2.745, 0.085,   12.0, 1.565},
{"B_3",   0.838, 109.47, 4.083, 0.180, 12.052, 1.755},
{"B_2",   0.828,  120.0, 4.083, 0.180, 12.052, 1.755},
{"C_3",   0.757, 109.47, 3.851, 0.105,  12.73, 1.912},
{"C_R",   0.729,  120.0, 3.851, 0.105,  12.73, 1.912},
{"C_2",   0.732,  120.0, 3.851, 0.105,  12.73, 1.912},
{"C_1",   0.706,  180.0, 3.851, 0.105,  12.73, 1.912},
{"N_3",   0.700,  106.7, 3.660, 0.069, 13.407, 2.544},
{"N_R",   0.699,  120.0, 3.660, 0.069, 13.407, 2.544},
{"N_2",   0.685,  111.2, 3.660, 0.069, 13.407, 2.544},
{"N_1",   0.656,  180.0, 3.660, 0.069, 13.407, 2.544},
{"O_3",   0.658, 104.51, 3.500, 0.060, 14.085, 2.300},
{"O_3_z", 0.528,  146.0, 3.500, 0.060, 14.085, 2.300},
{"O_R",   0.680,  110.0, 3.500, 0.060, 14.085, 2.300},
{"O_2",   0.634,  120.0, 3.500, 0.060, 14.085, 2.300},
{"O_1",   0.639,  180.0, 3.500, 0.060, 14.085, 2.300},
{"F_",    0.668,  180.0, 3.364, 0.050, 14.762, 1.735},
{"Ne4+4", 0.920,   90.0, 3.243, 0.042, 15.440, 0.194},
{"Na",    1.539,  180.0, 2.983, 0.030,   12.0, 1.081},
{"Mg3+2", 1.421, 109.47, 3.021, 0.111,   12.0, 1.787},
{"Al3",   1.244, 109.47, 4.499, 0.505, 11.278, 1.792},
{"Si3",   1.117, 109.47, 4.295, 0.402, 12.175, 2.323},
{"P_3+3", 1.101,   93.8, 4.147, 0.305, 13.072, 2.863},
{"P_3+5", 1.056, 109.47, 4.147, 0.305, 13.072, 2.863},
//{"P_3+q", 1.056, 109.47, 4.147, 0.305, 13.072, 2.863}, duplicate
{"S_3+2", 1.064,   92.1, 4.035, 0.274, 13.969, 2.703},
{"S_3+4", 1.049, 103.20, 4.035, 0.274, 13.969, 2.703},
{"S_3+6", 1.027, 109.47, 4.035, 0.274, 13.969, 2.703},
{"S_R",   1.077,   92.2, 4.035, 0.274, 13.969, 2.703},
{"S_2",   0.854,  120.0, 4.035, 0.274, 13.969, 2.703},
{"Cl",    1.044,  180.0, 3.947, 0.227, 14.866, 2.348},
{"Ar4+4", 1.032,   90.0, 3.868, 0.185, 15.763, 0.300},
{"K_",    1.953,  180.0, 3.812, 0.035,   12.0, 1.165},
{"Ca6+2", 1.761,   90.0, 3.399, 0.238,   12.0, 2.141},
{"Sc3+3", 1.513, 109.47, 3.295, 0.019,   12.0, 2.592},
{"Ti3+4", 1.412, 109.47, 3.175, 0.017,   12.0, 2.659},
{"Ti6+4", 1.412,   90.0, 3.175, 0.017,   12.0, 2.659},
{"V_3+5", 1.402, 109.47, 3.144, 0.016,   12.0, 2.679},
{"Cr6+3", 1.345,   90.0, 3.023, 0.015,   12.0, 2.463},
{"Mn6+2", 1.382,   90.0, 2.961, 0.013,   12.0,  2.43},
{"Fe3+2", 1.270, 109.47, 2.912, 0.013,   12.0,  2.43},
{"Fe6+2", 1.335,   90.0, 2.912, 0.013,   12.0,  2.43},
{"Co6+3", 1.241,   90.0, 2.872, 0.014,   12.0,  2.43},
{"Ni4+2", 1.164,   90.0, 2.834, 0.015,   12.0,  2.43},
{"Cu3+1", 1.302, 109.47, 3.495, 0.005,   12.0, 1.756},
{"Zn3+2", 1.193, 109.47, 2.763, 0.124,   12.0, 1.308},
{"Ga3+3", 1.260, 109.47, 4.383, 0.415,   11.0, 1.821},
{"Ge3",   1.197, 109.47, 4.280, 0.379,   12.0, 2.789},
{"As3+3", 1.211,   92.1, 4.230, 0.309,   13.0, 2.864},
{"Se3+2", 1.190,   90.6, 4.205, 0.291,   14.0, 2.764},
{"Br",    1.192,  180.0, 4.189, 0.251,   15.0, 2.519},
{"Kr4+4", 1.147,   90.0, 4.141, 0.220,   16.0, 0.452},
{"Rb",    2.260,  180.0, 4.114,  0.04,   12.0, 1.592},
{"Sr6+2", 2.052,   90.0, 3.641, 0.235,   12.0, 2.449},
{"Y_3+3", 1.698, 109.47, 3.345, 0.072,   12.0, 3.257},
{"Zr3+4", 1.564, 109.47, 3.124, 0.069,   12.0, 3.667},
{"Nb3+5", 1.473, 109.47, 3.165, 0.059,   12.0, 3.618},
{"Mo6+6", 1.467,   90.0, 3.052, 0.056,   12.0,  3.40},
{"Mo3+6", 1.484, 109.47, 3.052, 0.056,   12.0,  3.40},
{"Tc6+5", 1.322,   90.0, 2.998, 0.048,   12.0,  3.40},
{"Ru6+2", 1.478,   90.0, 2.963, 0.056,   12.0,  3.40},
{"Rh6+3", 1.332,   90.0, 2.929, 0.053,   12.0, 3.508},
{"Pd4+2", 1.338,   90.0, 2.899, 0.048,   12.0,  3.21},
{"Ag1+1", 1.386,  180.0, 3.148, 0.036,   12.0, 1.956},
{"Cd3+2", 1.403, 109.47, 2.848, 0.228,   12.0,  1.65},
{"In3+3", 1.459, 109.47, 4.463, 0.599,   11.0,  2.07},
{"Sn3",   1.398, 109.47, 4.392, 0.567,   12.0, 2.961},
{"Sb3+3", 1.407,   91.6, 4.420, 0.449,   13.0, 2.704},
{"Te3+2", 1.386,  90.25, 4.470, 0.398,   14.0, 2.882},
{"I_",    1.382,  180.0,  4.50, 0.339,   15.0,  2.65},
{"Xe4+4", 1.267,   90.0, 4.404, 0.332,   12.0, 0.556},
{"Cs",    2.570,  180.0, 4.517, 0.045,   12.0, 1.573},
{"Ba6+2", 2.277,   90.0, 3.703, 0.364,   12.0, 2.727},
{"La3+3", 1.943, 109.47, 3.522, 0.017,   12.0,  3.30},
{"Ce6+3", 1.841,   90.0, 3.556, 0.013,   12.0,  3.30},
{"Pr6+3", 1.823,   90.0, 3.606, 0.010,   12.0,  3.30},
{"Nd6+3", 1.816,   90.0, 3.575, 0.010,   12.0,  3.30},
{"Pm6+3", 1.801,   90.0, 3.547, 0.009,   12.0,  3.30},
{"Sm6+3", 1.780,   90.0, 3.520, 0.008,   12.0,  3.30},
{"Eu6+3", 1.771,   90.0, 3.493, 0.008,   12.0,  3.30},
{"Gd6+3", 1.735,   90.0, 3.368, 0.009,   12.0,  3.30},
{"Tb6+3", 1.732,   90.0, 3.451, 0.007,   12.0,  3.30},
{"Dy6+3", 1.710,   90.0, 3.428, 0.007,   12.0,  3.30},
{"Ho6+3", 1.696,   90.0, 3.409, 0.007,   12.0, 3.416},
{"Er6+3", 1.673,   90.0, 3.391, 0.007,   12.0,  3.30},
{"Tm6+3", 1.660,   90.0, 3.374, 0.006,   12.0,  3.30},
{"Yb6+3", 1.637,   90.0, 3.355, 0.228,   12.0, 2.618},
{"Lu6+3", 1.671,   90.0, 3.640, 0.041,   12.0, 3.271},
{"Hf3+4", 1.611, 109.47, 3.141, 0.072,   12.0, 3.921},
{"Ta3+5", 1.511, 109.47, 3.170, 0.081,   12.0, 4.075},
{"W_6+6", 1.392,   90.0, 3.069, 0.067,   12.0,  3.70},
{"W_3+4", 1.526, 109.47, 3.069, 0.067,   12.0,  3.70},
{"W_3+6", 1.380, 109.47, 3.069, 0.067,   12.0,  3.70},
{"Re6+5", 1.372,   90.0, 2.954, 0.066,   12.0,  3.70},
{"Re3+7", 1.314, 109.47, 2.954, 0.066,   12.0,  3.70},
{"Os6+6", 1.372,   90.0, 3.120, 0.037,   12.0,  3.70},
{"Ir6+3", 1.371,   90.0, 2.840, 0.073,   12.0, 3.731},
{"Pt4+2", 1.364,   90.0, 2.754, 0.080,   12.0, 3.382},
{"Au4+3", 1.262,   90.0, 3.293, 0.039,   12.0, 2.625},
{"Hgl+2", 1.340,  180.0, 2.705, 0.385,   12.0,  1.75},
{"T13+3", 1.518,  120.0, 4.347, 0.680,   11.0, 2.068},
{"Pb3",   1.459, 109.47, 4.297, 0.663,   12.0, 2.846},
{"Bi3+3", 1.512,   90.0, 4.370, 0.518,   13.0, 2.470},
{"Po3+2",  1.50,   90.0, 4.709, 0.325,   14.0,  2.33},
{"At",    1.545,  180.0, 4.750, 0.284,   15.0,  2.24},
{"Rn4+4", 1.420,   90.0, 4.765, 0.248,   16.0, 0.583},
{"Fr",    2.880,  180.0,  4.90, 0.050,   12.0, 1.847},
{"Ra6+2", 2.512,   90.0, 3.677, 0.404,   12.0,  2.92},
{"Ac6+3", 1.983,   90.0, 3.478, 0.033,   12.0,  3.90},
{"Th6+4", 1.721,   90.0, 3.396, 0.026,   12.0, 4.202},
{"Pa6+4", 1.711,   90.0, 3.424, 0.022,   12.0,  3.90},
{"U_6+4", 1.684,   90.0, 3.395, 0.022,   12.0,  3.90},
{"Np6+4", 1.666,   90.0, 3.424, 0.019,   12.0,  3.90},
{"Pu6+4", 1.657,   90.0, 3.424, 0.016,   12.0,  3.90},
{"Am6+4", 1.660,   90.0, 3.381, 0.014,   12.0,  3.90},
{"Cm6+3", 1.801,   90.0, 3.326, 0.013,   12.0,  3.90},
{"Bk6+3", 1.761,   90.0, 3.339, 0.013,   12.0,  3.90},
{"Cf6+3", 1.750,   90.0, 3.313, 0.013,   12.0,  3.90},
{"Es6+3", 1.724,   90.0, 3.299, 0.012,   12.0,  3.90},
{"Fm6+3", 1.712,   90.0, 3.286, 0.012,   12.0,  3.90},
{"Md6+3", 1.689,   90.0, 3.274, 0.011,   12.0,  3.90},
{"No6+3", 1.679,   90.0, 3.248, 0.011,   12.0,  3.90},
{"Lw6+3", 1.698,   90.0, 3.236, 0.011,   12.0,  3.90}};
