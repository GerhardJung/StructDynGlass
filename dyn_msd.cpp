
#include "dyn_msd.h"
#include "global.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_msd(){

    int flag = msd_flag;

    // loop over time
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL MSD " << t << std::endl; 
        reset_dyn(t);

        if (dyn_ranges[flag][0] > -0.001) {
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
                    dyn_loc[s*N+i] = dr;
                }
            }
            calc_bonds(dyn_loc,&dyn_ranges_time[flag][2*t]);
        }

        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) {
                for (int j=0; j<NI;j++) {
                    double dr = 0, dx;
                    for (int d=0; d<dim;d++) {
                        dx = xyz_data[i+s*N][d+dim*NT*j] - xyz_data[i+s*N][d+t*dim+dim*NT*j];
                        apply_pbc(dx);
                        dr += dx*dx;
                    }
                    add_histogram_avg(s,i,dyn_ranges_time[flag][2*t],dyn_ranges_time[flag][2*t+1],dr);
                }
            }
        }

        // the main evaluation for the isoconfigurational ensemble
        eval_isoconf(t, flag);

    }

    // write results
    print_isoconf(flag);

}