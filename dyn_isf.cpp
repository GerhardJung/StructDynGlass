
#include "dyn_isf.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_isf(){

    int flag = isf_flag;

    double * save_bb_traj = new double [NS*NT*NI*N];

    int d=0;
    // loop over time
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL ISF " << t << std::endl; 
        reset_dyn(t);

        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) {
                for (int j=0; j<NI;j++) {
                    double dx = xyz_data[i+s*N][d+NDim*NT*j] - xyz_data[i+s*N][d+t*NDim+NDim*NT*j];
                    apply_pbc(dx);
                    double C_loc= cos(qisf*dx);

                    add_histogram_avg(s,i,j,dyn_ranges[flag],C_loc);
                    save_bb_traj[s*NT*NI*N+t*N*NI+j*N+i] =C_loc;
                }
            }
        }

        // the main evaluation for the isoconfigurational ensemble
        eval_isoconf(t,flag);
    }
    
    // write results
    print_traj(save_bb_traj, flag);

    // calculate and write R file
    print_R(save_bb_traj, flag);

    // write results
    print_isoconf(isf_flag);
}

