
#include "dyn_rearrange_patinet.h"
#include "struct_soft_modes.h"
#include "global.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"
#include <stdlib.h>     /* srand, rand */

void eval_rp(){

    int flag = rp_flag;

    // loop over time
    save_pat = new double [3*NS*NI*N];
    save_pat_traj = new double [NT*NI*N];
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL REARRANGE PATINET " << t << std::endl; 

        if (dyn_ranges[flag][0] > -10000.0) {
            dyn_ranges_time[flag][3*t] = dyn_ranges[flag][0];
            dyn_ranges_time[flag][3*t+1] = dyn_ranges[flag][1];
            dyn_ranges_time[flag][3*t+2] = dyn_ranges[flag][2];
        } else { // determine binning from maximal and minimal values
            // TODO
        }
        if (dyn_ranges[flag+1][0] > -10000.0) {
            dyn_ranges_time[flag+1][3*t] = dyn_ranges[flag+1][0];
            dyn_ranges_time[flag+1][3*t+1] = dyn_ranges[flag+1][1];
            dyn_ranges_time[flag+1][3*t+2] = dyn_ranges[flag+1][2];
        } else { // determine binning from maximal and minimal values
            // TODO
        }
        if (dyn_ranges[flag+2][0] > -10000.0) {
            dyn_ranges_time[flag+2][3*t] = dyn_ranges[flag+2][0];
            dyn_ranges_time[flag+2][3*t+1] = dyn_ranges[flag+2][1];
            dyn_ranges_time[flag+2][3*t+2] = dyn_ranges[flag+2][2];
        } else { // determine binning from maximal and minimal values
            // TODO
        }

        // first evaluate quantities
        double uth[N*dim] = {0};
        double flin[N*dim] = {0};
          srand (time(NULL));
        for (int s=0; s<NS;s++) { // loop over structures

            if (dyn_rearrange_mode==0) {
                double result_loc[dim*dim], dx[dim];
                for (int i=0; i<N;i++) { // loop over particles
                    for (int i2=0; i2<N;i2++) { // loop over particles
                        calc_2Depot(i+s*N,i2+s*N,0,0,result_loc,dx,hessian[i*N+i2]);
                    }
                }
            }


            for (int j=0; j<NI;j++) {

                //std::cout << j << " START HESSIAN" << std::endl;

                if (dyn_rearrange_mode==1) {
                    double result_loc[dim*dim], dx[dim];
                    for (int i=0; i<N;i++) { // loop over particles
                        for (int i2=0; i2<N;i2++) { // loop over particles
                            calc_2Depot(i+s*N,i2+s*N,t,j,result_loc,dx,hessian[i*N+i2]);
                        }
                    }
                }

                //std::cout << j << " FIN HESSIAN" << std::endl;


                for (int i=0; i<N;i++) {
                    for (int d=0; d<dim;d++) {
                        uth[d+i*dim] = xyz_inherent_data[i+s*N][d+dim*NT*j] - xyz_inherent_data[i+s*N][d+t*dim+dim*NT*j];
                        //std::cout << xyz_inherent_data[i+s*N][d+dim*NT*j] << " " << xyz_inherent_data[i+s*N][d+t*dim+dim*NT*j] << std::endl;
                        apply_pbc(uth[d+i*dim]);
                        //uth[d+i*dim] = ((double) rand() / (RAND_MAX)) + -0.5;
                    }
                    //std::cout << "diff: " << uth[0+i*dim] << " " << uth[1+i*dim] << std::endl;
                }
                // multiply with hessian
                for (int i=0; i<N*dim;i++) flin[i] = 0.0;
                for (int i=0; i<N;i++) {
                    for (int di=0; di<dim;di++) {
                        for (int i2=0; i2<N;i2++) {
                            for (int di2=0; di2<dim;di2++) {
                                flin[di+i*dim] += hessian[i*N+i2][di2+di*dim]*uth[di2+i2*dim];
                                //std::cout << "hessian " << hessian[s*N*N+i*N+i2][di2+di*dim] << std::endl;
                            }
                        }
                    }
                    double uth2 = 0.0;
                    for (int di=0; di<dim;di++) uth2 += uth[di+i*dim]*uth[di+i*dim];
                    uth2 = sqrt(uth2);
                    if (uth2 < 1e-6) uth2 = 1e-6;
                    double log10uth2 = log10(uth2);
                    save_pat[s*N*NI+j*N+i] = log10uth2;

                    double fres = 0.0;
                    for (int di=0; di<dim;di++) {
                        fres += flin[di+i*dim]*flin[di+i*dim];
                        //std::cout << "flin " << flin[di+i*dim] << std::endl;
                    }
                    //std::cout << fres << std::endl;
                    fres = sqrt(fres);
                    if (fres < 1e-6) fres = 1e-6;
                    double passive = 1.0;
                    if (fres > dyn_rearrange_threshold) passive = 0.0;
                    double log10fres = log10(fres);
                    save_pat[NS*NI*N+s*N*NI+j*N+i] = log10fres;
                    //save_pat[NS*NI*N+s*N*NI+j*N+i] = flin[i*dim];
                    save_pat[2*NS*NI*N+s*N*NI+j*N+i] = passive;
                    if (s==0) save_pat_traj[t*N*NI+j*N+i] = log10fres;
                }

                //std::cout << j << " FIN REST" << std::endl;
            }
        }

        //std::cout << save_pat[0*NS*NI*N+1*N*NI+10*N+10] << " " << save_pat[1*NS*NI*N+1*N*NI+10*N+10] << " " << save_pat[2*NS*NI*N+1*N*NI+10*N+10] << "\n";

        // then save histograms for all three dynamical observables
        for (int k=0; k<3; k++) {
            reset_dyn(t);
            //if (k==1) std::cout << dyn_ranges_time[flag+k][2*t]
            for (int s=0; s<NS;s++) { // loop over structures
                for (int j=0; j<NI;j++) {
                    for (int i=0; i<N;i++) {
                        add_histogram_avg(s,i,j,&dyn_ranges_time[flag+k][3*t],save_pat[k*NS*NI*N+s*N*NI+j*N+i]);
                        //if (t==2) std::cout << fres << " " << dyn_ranges_time[flag][2*t] << " "  << dyn_ranges_time[flag][2*t+1] << std::endl;
                    }
                }
            }
            // the main evaluation for the isoconfigurational ensemble
            eval_isoconf(t, flag+k);
        }

    }

    // write results
    print_traj(save_pat_traj,flag);
    print_isoconf(flag);
    print_isoconf(flag+1);
    print_isoconf(flag+2);

}

