
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
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL REARRANGE PATINET " << t << std::endl; 

        if (dyn_ranges[flag][0] > -10000.0) {
            dyn_ranges_time[flag][2*t] = dyn_ranges[flag][0];
            dyn_ranges_time[flag][2*t+1] = dyn_ranges[flag][1];
        } else { // determine binning from maximal and minimal values
            // TODO
        }
        if (dyn_ranges[flag+1][0] > -10000.0) {
            dyn_ranges_time[flag+1][2*t] = dyn_ranges[flag+1][0];
            dyn_ranges_time[flag+1][2*t+1] = dyn_ranges[flag+1][1];
        } else { // determine binning from maximal and minimal values
            // TODO
        }
        if (dyn_ranges[flag+2][0] > -10000.0) {
            dyn_ranges_time[flag+2][2*t] = dyn_ranges[flag+2][0];
            dyn_ranges_time[flag+2][2*t+1] = dyn_ranges[flag+2][1];
        } else { // determine binning from maximal and minimal values
            // TODO
        }

        // first evaluate quantities
        double uth[N*dim] = {0};
        double flin[N*dim] = {0};
          srand (time(NULL));
        for (int s=0; s<NS;s++) { // loop over structures
            for (int j=0; j<NI;j++) {

                // calculate hessian
                if (struct_soft_modes_flag==-1) {
                    for (int i=0; i<N;i++) { // loop over particles
                        for (int i2=0; i2<N;i2++) { // loop over particles
                            calc_2Depot(i+s*N,i2+s*N,t,j,hessian[s*N*N+i*N+i2]);
                        }
                    }
                }


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
                                flin[di+i*dim] += hessian[s*N*N+i*N+i2][di2+di*dim]*uth[di2+i2*dim];
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
                    for (int di=0; di<dim;di++) fres += flin[di+i*dim]*flin[di+i*dim];
                    fres = sqrt(fres);
                    if (fres < 1e-6) fres = 1e-6;
                    double passive = 1.0;
                    if (fres > dyn_rearrange_threshold) passive = 0.0;
                    double log10fres = log10(fres);
                    save_pat[NS*NI*N+s*N*NI+j*N+i] = log10fres;
                    save_pat[2*NS*NI*N+s*N*NI+j*N+i] = passive;

                }
            }
        }

        //std::cout << save_pat[0*NS*NI*N+1*N*NI+10*N+10] << " " << save_pat[1*NS*NI*N+1*N*NI+10*N+10] << " " << save_pat[2*NS*NI*N+1*N*NI+10*N+10] << "\n";

        // then save histograms for all three dynamical observables
        for (int k=0; k<3; k++) {
            reset_dyn(t);
            for (int s=0; s<NS;s++) { // loop over structures
                for (int j=0; j<NI;j++) {
                    for (int i=0; i<N;i++) {
                        add_histogram_avg(s,i,dyn_ranges_time[flag+k][2*t],dyn_ranges_time[flag+k][2*t+1],save_pat[k*NS*NI*N+s*N*NI+j*N+i]);
                        //if (t==2) std::cout << fres << " " << dyn_ranges_time[flag][2*t] << " "  << dyn_ranges_time[flag][2*t+1] << std::endl;
                    }
                }
            }
            // the main evaluation for the isoconfigurational ensemble
            eval_isoconf(t, flag+k);
        }

    }

    // write results
    //print_traj(save_pat);
    print_isoconf(flag);
    print_isoconf(flag+1);
    print_isoconf(flag+2);

}

