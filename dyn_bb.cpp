
#include "dyn_bb.h"
#include "defs.h"
#include "pbc.h"
#include "struct_base.h"
#include "eval_isoconf.h"
#include "eval_struct.h"

void eval_bb(){

    int flag = bb_flag;

    double * save_bb_traj = new double [NS*NT*NI*N];
    int ** neighbors = imatrix(0,NS*N-1,0,N_NEIGH_MAX-1);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            for (int d=0; d<N_NEIGH_MAX;d++) {
                neighbors[i+s*N][d] = -1;
            }
        }
    }

    // calc original neighbors
    findneighbors(rcuti2, neighbors);

    // loop over time
    int n0;
    int nt;
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL BB " << t << std::endl; 
        reset_dyn(t);

        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) {
                for (int j=0; j<NI;j++) {
                    checkneighbors(s, i,j, t, n0, nt, neighbors);
                    // add to probability distribution and averages
                    double C_loc= nt/((double) n0);
                    add_histogram_avg(s,i,j,dyn_ranges[flag],C_loc);
                    //if (t==1 && C_loc <1) std::cout << s << " "<< i << " " << C_loc << std::endl;
                    save_bb_traj[s*NT*NI*N+t*N*NI+j*N+i] =C_loc;
                }
            }
        }

        // the main evaluation for the isoconfigurational ensemble
        eval_isoconf(t,flag);
    }

    // calculate rearranging time scale
    eval_timescale(flag, 0.5);


    // just S and G
    if (boxL > 100.0) {
        double tmp[N*NS];
        for (int s=0; s<NS; s++) {
            for (int i=0; i<N; i++) {
                tmp[i] = 1.0;
            }
        }
        std::string tmps = "BASE";
        int first =1;
        eval_struct(tmp,tmps,first);

        // calc S4, G4
        for (int t=15; t<NT; t++) {
            first = 0;

            for (int i=0; i<N*NS; i++) {
                tmp[i] = 1.0-dyn_avg_save[i][flag*(NT+1)+t];
            }
            tmps = "BB"+std::to_string(t);
            std::cout << tmps << " " << tmp[2] << std::endl;
            eval_struct(tmp,tmps,first);
            for (int j=0; j<NI;j++) {
                for (int s=0; s<NS; s++) {
                    for (int i=0; i<N; i++) {
                        tmp[i+s*N] = 1.0-save_bb_traj[s*NT*NI*N+t*N*NI+j*N+i];
                    }
                }
                tmps = "BB"+std::to_string(t)+"Traj"+std::to_string(j);
                eval_struct(tmp,tmps,first);
            }

        }
    }

        


    // write results
    print_traj(save_bb_traj, flag);

    // write results
    print_isoconf(flag);

    free_imatrix(neighbors,0,NS*N-1,0,N_NEIGH_MAX-1);
}



// Help functions
void findneighbors(double rcut2, int ** neighbors) {
    double dr, dx;
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            int ncount = 0;
            for (int j=0; j<N;j++) {
                dr = 0;
                for (int d=0; d<NDim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                }

                double sigma = determine_sigma(i,j);

                if (dr < rcut2*sigma*sigma && i != j ) {
                    neighbors[i+s*N][ncount] = j;
                    ncount ++;
                }
            }
            //if (s==0 && i == 264) std::cout << rcut2 << std::endl;
        }
    }
}

void checkneighbors(int s, int i, int j, int t, int &n0, int &nt, int ** neighbors) {
    double dr, dx;
    n0 = 0;
    nt = 0;
    while (neighbors[i+s*N][n0] != -1 ) {
        //if(n0>5) std::cout << neighbors[i+s*N][n0] << std::endl;
        dr = 0;
        for (int d=0; d<NDim;d++) {
            dx = xyz_data[i+s*N][d+t*NDim+NDim*NT*j] - xyz_data[neighbors[i+s*N][n0]+s*N][d+t*NDim+NDim*NT*j];
            apply_pbc(dx);
            dr += dx*dx;
        }

        double sigma = determine_sigma(i,j);

        if (dr < rcuto2*sigma*sigma) {
            nt ++;
        }
        n0 ++;
    }
}