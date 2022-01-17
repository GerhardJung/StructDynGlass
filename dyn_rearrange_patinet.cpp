
#include "dyn_rearrange_patinet.h"
#include "struct_soft_modes.h"
#include "global.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_rp(){

    int flag = rp_flag;

    if (struct_soft_modes_flag==-1) {
        // calculate hessian
        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) { // loop over particles
                for (int j=0; j<N;j++) { // loop over particles
                    calc_2Depot(i+s*N,j+s*N,hessian[s*N*N+i*N+j]);
                }
            }
        }
    }

    // loop over time
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL REARRANGE PATINET " << t << std::endl; 
        reset_dyn(t);

        /*if (dyn_ranges[flag][0] > -0.001) {
            dyn_ranges_time[flag][2*t] = dyn_ranges[flag][0];
            dyn_ranges_time[flag][2*t+1] = dyn_ranges[flag][1];
        } else { // determine binning from maximal and minimal values
            double dyn_loc[N*NS] = {0};
            for (int s=0; s<NS;s++) { // loop over structures
                for (int i=0; i<N;i++) {
                    double dr = 0, dx;
                    for (int d=0; d<dim;d++) {
                        dx = xyz_data[i+s*N][d] - xyz_data[i+s*N][d+t*dim];
                        apply_pbc(dx);
                        dr += dx*dx;
                    }
                    dyn_loc[s*N+i] = exp(-dr*dr*exp_scale4i);
                }
            }
            calc_bonds(dyn_loc,&dyn_ranges_time[flag][2*t]);
        }*/

        double uth[N*dim] = {0};
        double flin[N*dim] = {0};

        for (int s=0; s<NS;s++) { // loop over structures
            for (int j=0; j<NI;j++) {
                for (int i=0; i<N;i++) {
                    for (int d=0; d<dim;d++) {
                        uth[d+i*dim] = xyz_inherent_data[i+s*N][d+dim*NT*j] - xyz_inherent_data[i+s*N][d+t*dim+dim*NT*j];
                        apply_pbc(uth[d+i*dim]);
                    }
                }
                // multiply with hessian
                for (int i=0; i<N*dim;i++) flin[i] = 0.0;
                for (int i=0; i<N;i++) {
                    for (int di=0; di<dim;di++) {
                        for (int j=0; j<N;j++) {
                            for (int dj=0; dj<dim;dj++) {
                                flin[di+i*dim] += hessian[s*N*N+i*N+j][dj+di*dim]*uth[dj+j*dim];
                            }
                        }
                    }
                    double fres = 0.0;
                    for (int di=0; di<dim;di++) fres += flin[di+i*dim]*flin[di+i*dim];
                    fres = sqrt(fres);
                    add_histogram_avg(s,i,dyn_ranges_time[flag][2*t],dyn_ranges_time[flag][2*t+1],fres);
                }

            }
        }

        // the main evaluation for the isoconfigurational ensemble
        eval_isoconf(t, flag);

    }

    // write results
    print_isoconf(flag);


}

