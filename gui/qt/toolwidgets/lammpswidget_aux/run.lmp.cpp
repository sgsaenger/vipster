#include "vipster/global.h"

#include "run.lmp.h"
#include "toolwidgets/lammpswidget_aux/fix_vipster.lmp.h"

#include <unistd.h>

#include <chrono>
#include <thread>
#include <iostream>
#include <QApplication>
#include <fmt/format.h>

#include "lammps.h"
#include "comm.h"
#include "input.h"
#include "modify.h"
#include "exceptions.h"
#include "error.h"

using namespace std::chrono_literals;
using namespace LAMMPS_NS;
using namespace Vipster;

static std::atomic<bool> running{false};

std::pair<int, std::string> Lammps::runMaster(std::string dir, runParams params, Molecule *mol)
{
#ifdef USE_MPI
    if(params.MPI <= 1){
#endif
        // sequential execution on main vipster process
        try{
            run(dir, params, MPI_COMM_NULL, mol);
        }catch(std::exception &e){
            return {1, e.what()};
        }
        return {};
#ifdef USE_MPI
    }else{
        running = true;
        // parallel execution via spawned child processes
        char *argv[2];
        argv[0] = "lammps_mpi_slave";
        argv[1] = nullptr;
        MPI_Info info;
        MPI_Info_create(&info);
        MPI_Info_set(info, "map_by", "core:OVERSUBSCRIBE");
        MPI_Comm intercomm;
        int spawned = MPI_Comm_spawn(
            QApplication::arguments()[0].toUtf8(),
            argv,
            params.MPI,
            info,
            0,
            MPI_COMM_SELF,
            &intercomm,
            MPI_ERRCODES_IGNORE);
        MPI_Info_free(&info);
        if(spawned){
            // failed to create child processes
            throw Vipster::Error{"Failed to launch LAMMPS processes."};
        }
        // send working directory
        unsigned long len_dir = dir.size()+1; // also submit terminating null char
        MPI_Bcast(&len_dir, 1, MPI_UNSIGNED_LONG, MPI_ROOT, intercomm);
        MPI_Bcast(&dir[0], len_dir, MPI_CHAR, MPI_ROOT, intercomm);
        // send parameters
        MPI_Bcast(&params, sizeof(runParams), MPI_BYTE, MPI_ROOT, intercomm);
        // wait for lammps to finish
        std::pair<int, std::string> res;
        while(running){
            int flag{0};
            MPI_Status status{};
            MPI_Iprobe(MPI_ANY_SOURCE, 0, intercomm, &flag, &status);
            if(flag){
                MPI_Recv(&res.first, 1, MPI_INT, status.MPI_SOURCE, 0, intercomm, &status);
                if(res.first){
                    // lammps failed, received error message
                    long msg{};
                    MPI_Recv(&msg, 1, MPI_LONG, status.MPI_SOURCE, 1, intercomm, &status);
                    // lammps failed, receive error message
                    res.second.resize(msg, ' ');
                    MPI_Recv(&res.second[0], msg, MPI_CHAR, status.MPI_SOURCE, 2, intercomm, &status);
                }else{
                    // lammps succeeded, return without message
                }
                running = false;
            }
        }
        // TODO: handle copying Lammps data to GUI
        MPI_Comm_disconnect(&intercomm);
        return res;
    }
}

void Lammps::runSlave()
{
    running = true;
    int argc{0};
    char **argv{nullptr};
    int level{0};
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
    MPI_Comm intercomm;
    MPI_Comm_get_parent(&intercomm);
    int me{0};
    int size{0};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    // get working directory
    std::string dir;
    unsigned long len_dir;
    MPI_Bcast(&len_dir, 1, MPI_UNSIGNED_LONG, 0, intercomm);
    char buf[len_dir];
    MPI_Bcast(buf, len_dir, MPI_CHAR, 0, intercomm);
    dir = buf;
    // get run parameters
    runParams params;
    MPI_Bcast(&params, sizeof(runParams), MPI_BYTE, 0, intercomm);
    // handle abort events via a seperate thread
    std::thread error_handler{[&](){
        int flag{0};
        MPI_Status status{};
        while(running){
            MPI_Iprobe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &flag, &status);
            if(flag){
                MPI_Recv(&flag, 1, MPI_INT, status.MPI_SOURCE, 3, MPI_COMM_WORLD, &status);
                MPI_Comm_disconnect(&intercomm);
                MPI_Finalize();
                exit(0);
            }
            std::this_thread::sleep_for(1s);
        }
    }};
    // actual run
    int success{0};
    int error{1};
    int fatal{-1};
    try{
        run(dir, params, intercomm);
        if(me == 0){
            MPI_Send(&success, 1, MPI_INT, 0, 0, intercomm);
        }
    }catch(LAMMPSAbortException &e){
        // single-process failure, abort slave processes
        long msgSize = e.message.size();
        MPI_Send(&fatal, 1, MPI_INT, 0, 0, intercomm);
        MPI_Send(&msgSize, 1, MPI_LONG, 0, 1, intercomm);
        MPI_Send(&e.message[0], e.message.size(), MPI_CHAR, 0, 2, intercomm);
        int abort{0};
        for(int i=0; i<size; ++i){
            if(i == me) continue;
            MPI_Send(&abort, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
        }
    }catch(LAMMPSException &e){
        // synchronous failure, finish normally
        if(me == 0){
            long msgSize = e.message.size();
            MPI_Send(&error, 1, MPI_INT, 0, 0, intercomm);
            MPI_Send(&msgSize, 1, MPI_LONG, 0, 1, intercomm);
            MPI_Send(&e.message[0], e.message.size(), MPI_CHAR, 0, 2, intercomm);
        }
    }catch(...){
        // unknown failure, abort slave processes
        const char* msg = "Unknown error in LAMMPS calculation.";
        long msgSize = 36;
        MPI_Send(&fatal, 1, MPI_INT, 0, 0, intercomm);
        MPI_Send(&msgSize, 1, MPI_LONG, 0, 1, intercomm);
        MPI_Send(msg, msgSize, MPI_CHAR, 0, 2, intercomm);
        int abort{0};
        for(int i=0; i<size; ++i){
            if(i == me) continue;
            MPI_Send(&abort, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
        }
    }
    running = false;
    error_handler.join();
    MPI_Comm_disconnect(&intercomm);
    MPI_Finalize();
#endif
}

void Lammps::run(std::string dir, runParams params, MPI_Comm intercomm, Molecule *mol)
{
    if(!(intercomm == MPI_COMM_NULL) ^ (mol == nullptr)){
        throw Vipster::Error{"Only use MPI or in-memory access to Stepdata."};
    }
    auto log = dir+"/log.lammps";
    char* lmparg[5]{
        nullptr,
        "-screen", "none",
        "-log", &log[0]
    };
    LAMMPS lmp{5, lmparg, MPI_COMM_WORLD};
    // register custom fix
    (*lmp.modify->fix_map)["vipster"] = &mkFixVipster;
    // setup lammps
    lmp.input->file((dir+"/input").c_str());
    // initialize fix vipster
    if(!lmp.modify || !lmp.modify->fix){
        throw Vipster::Error{"Lammps could not be initizialized succesfully."};
    }
    auto fix_vipster = dynamic_cast<FixVipster*>(lmp.modify->fix[lmp.modify->nfix-1]);
    if(!fix_vipster){
        throw Vipster::Error{"Error on registering callback fix."};
    }
    // register callback
    fix_vipster->mol = mol;
    fix_vipster->intercomm = intercomm;
    if(params.mode == runParams::Mode::MD){
        lmp.input->one(fmt::format("run {}", params.nstep).c_str());
    }else if(params.mode == runParams::Mode::Min){
        lmp.input->one(fmt::format("minimize {} {} {} {}",
                                   params.etol,
                                   params.ftol,
                                   params.nstep,
                                   params.neval).c_str());
    }else{
        throw Vipster::Error{fmt::format("Invalid Lammps run mode {} on proc {}.", params.mode, lmp.comm->me)};
    }
}
