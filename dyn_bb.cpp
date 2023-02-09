
#include "dyn_bb.h"
#include "defs.h"
#include "pbc.h"
#include "struct_base.h"
#include "eval_isoconf.h"
#include "eval_struct.h"

void eval_bb(){

    int flag = bb_flag;

    printf("385 (t=0) %f/%f, 385 (t=0.03) %f/%f,993 (t=0) %f/%f, 993 (t=0.03) %f/%f \n",xyz_data[385][0+0*NDim+NDim*NT*7],\
    xyz_data[385][1+0*NDim+NDim*NT*7],xyz_data[385][0+3*NDim+NDim*NT*7],xyz_data[385][1+3*NDim+NDim*NT*7],xyz_data[993][0+0*NDim+NDim*NT*7],\
    xyz_data[993][1+0*NDim+NDim*NT*7],xyz_data[993][0+3*NDim+NDim*NT*7],xyz_data[993][1+3*NDim+NDim*NT*7]);

    printf(" (t=0) %f, (t=0.03) %f\n",sqrt((xyz_data[385][0+0*NDim+NDim*NT*7]-xyz_data[993][0+0*NDim+NDim*NT*7])*(xyz_data[385][0+0*NDim+NDim*NT*7]-xyz_data[993][0+0*NDim+NDim*NT*7])+(xyz_data[385][1+0*NDim+NDim*NT*7]-xyz_data[993][1+0*NDim+NDim*NT*7])*(xyz_data[385][1+0*NDim+NDim*NT*7]-xyz_data[993][1+0*NDim+NDim*NT*7])),\
    sqrt((xyz_data[385][0+3*NDim+NDim*NT*7]-xyz_data[993][0+3*NDim+NDim*NT*7])*(xyz_data[385][0+3*NDim+NDim*NT*7]-xyz_data[993][0+3*NDim+NDim*NT*7])+(xyz_data[385][1+3*NDim+NDim*NT*7]-xyz_data[993][1+3*NDim+NDim*NT*7])*(xyz_data[385][1+3*NDim+NDim*NT*7]-xyz_data[993][1+3*NDim+NDim*NT*7]))
    );

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
                    //if (t==1 ) std::cout << "Test" << std::endl;
                    checkneighbors(s, i,j, t, n0, nt, neighbors, flag);
                    // add to probability distribution and averages
                    double C_loc= nt/((double) n0);
                    add_histogram_avg(s,i,j,dyn_ranges[flag],C_loc);
                    if (t==1 && C_loc <1) std::cout << s << " "<< i << " " << C_loc << std::endl;
                    save_bb_traj[s*NT*NI*N+t*N*NI+j*N+i] =C_loc;
                    //if (t==3 && C_loc < 1) printf("%d %d %f\n", i,j,C_loc); 
                    
                }
            }
        }

        // the main evaluation for the isoconfigurational ensemble
        eval_isoconf(t,flag);


#ifdef USE_RELATIVE
    for (int i=0; i<N*NS; i++) {
        for (int j=0; j<N_NEIGH_MAX;j++) {
            dyn_avg_save[i][(flag+j+1)*(NT+1)+t] /= (double) NI;
        }
    }
#endif
    }

    std::cout << "test1" << std::endl;
    // calculate rearranging time scale
    eval_timescale(flag, 0.5);
    // calc and save histograms for actual motion and isoconfigurational averages
    double hist_lower = hist_lower_time;
    double hist_upper = hist_upper_time;
    for (int i = 0; i < NS*N; i++) {
        int valint_dyn;
        double val = log10 (dyn_avg_save[i][flag*(NT+1)+NT]);
        if (val > hist_upper - EPS && val < hist_upper + EPS) valint_dyn = NHisto - 1;
        else valint_dyn = (val-hist_lower)/(hist_upper - hist_lower)* ((double)NHisto);

        if(valint_dyn >= 0 && valint_dyn < NHisto) dyn_hist_iso[NT+(NT+1)*flag][valint_dyn+NHisto*type_data[i]]+= 1.0;

    }

    // normalize histograms
    for (int l = 0; l < NHisto; l++) {
        for (int type=0; type<NTYPE; type++) {
            dyn_hist_iso[NT+(NT+1)*flag][l+NHisto*type] /= (double) NS*NPerType[type];
        }
    }

    // just S and G
    if (boxL > 100.0 && !noinherent) {
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
        for (int t=35; t<37; t++) {
            first = 0;

            for (int i=0; i<N*NS; i++) {
                tmp[i] = dyn_avg_save[i][flag*(NT+1)+t];
            }
            tmps = "BB"+std::to_string(t);
            
            eval_struct(tmp,tmps,first);
            /*for (int j=0; j<NI;j++) {
                for (int s=0; s<NS; s++) {
                    for (int i=0; i<N; i++) {
                        tmp[i+s*N] = save_bb_traj[s*NT*NI*N+t*N*NI+j*N+i];
                    }
                }
                tmps = "BB"+std::to_string(t)+"Traj"+std::to_string(j);
                eval_struct(tmp,tmps,first);
            }*/

        }
    }


    // write trajectories
    print_traj(save_bb_traj, flag);

    // write results
    print_isoconf(flag);

    // calculate and write R file
    print_R(save_bb_traj, flag);

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
                //double sigma = 1.0;

                if (dr < rcut2*sigma*sigma && i != j ) {
                    if (i==385) printf("%d %f %f\n",j,sqrt(dr),sigma);
                    neighbors[i+s*N][ncount] = j;
                    ncount ++;
                    if (ncount == N_NEIGH_MAX) printf("ERROR: Increase NMax for neighbor list! (BB)");
                }
            }
            //if (s==0 && i == 264) std::cout << rcut2 << std::endl;
        }
    }
}

void checkneighbors(int s, int i, int j, int t, int &n0, int &nt, int ** neighbors, int flag) {
    double dr, dx;
    n0 = 0;
    nt = 0;
    while (neighbors[i+s*N][n0] != -1 ) {
        //if(n0>5) std::cout << neighbors[i+s*N][n0] << std::endl;
        dr = 0.0;
        for (int d=0; d<NDim;d++) {
            dx = xyz_data[i+s*N][d+t*NDim+NDim*NT*j] - xyz_data[neighbors[i+s*N][n0]+s*N][d+t*NDim+NDim*NT*j];
            apply_pbc(dx);
            dr += dx*dx;
        }

        double sigma = determine_sigma(i,neighbors[i+s*N][n0]);
        //double sigma = 1.0;

        if (dr < rcuto2*sigma*sigma) {
            nt ++;
        }

        #ifdef USE_RELATIVE
            double dr0 = 0.0;
            for (int d=0; d<NDim;d++) {
                dx = xyz_data[i+s*N][d] - xyz_data[neighbors[i+s*N][n0]+s*N][d];
                apply_pbc(dx);
                dr0 += dx*dx;
            }
            dyn_avg_save[i][(flag+n0+1)*(NT+1)+t] +=  dr - dr0;
        #endif

        n0 ++;
    }
}