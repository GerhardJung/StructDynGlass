
#include "dyn_msd.h"
#include "global.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_msd(){

    int flag = msd_flag;

    // loop over time
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL MD " << t << std::endl; 

        // first calculate bonds (by just iterating over one trajectory)
        if (dyn_ranges[flag][0] > -10000.0) {
            dyn_ranges_time[flag][2*t] = dyn_ranges[flag][0];
            dyn_ranges_time[flag][2*t+1] = dyn_ranges[flag][1];
        } else { // determine binning from maximal and minimal values (just iterate over one simulation)
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
        if (dyn_ranges[flag+1][0] > -10000.0) {
            dyn_ranges_time[flag+1][2*t] = dyn_ranges[flag+1][0];
            dyn_ranges_time[flag+1][2*t+1] = dyn_ranges[flag+1][1];
        } else { // determine binning from maximal and minimal values (just iterate over one simulation)
            double dyn_loc[N*NS] = {0};
            for (int s=0; s<NS;s++) { // loop over structures
                for (int i=0; i<N;i++) {
                    double dr = 0, dx;
                    for (int d=0; d<dim;d++) {
                        dx = xyz_data[i+s*N][d] - xyz_data[i+s*N][d+t*dim];
                        apply_pbc(dx);
                        dr += dx*dx;
                    }
                    dyn_loc[s*N+i] = log(dr);
                }
            }
            calc_bonds(dyn_loc,&dyn_ranges_time[flag+1][2*t]);
        }

        // the main evaluation for the isoconfigurational ensemble (mean displacement)
        reset_dyn(t);
        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) {
                for (int j=0; j<NI;j++) {
                    double dr = 0, dx;
                    for (int d=0; d<dim;d++) {
                        dx = xyz_data[i+s*N][d+dim*NT*j] - xyz_data[i+s*N][d+t*dim+dim*NT*j];
                        apply_pbc(dx);
                        dr += dx*dx;
                    }
                    dr = sqrt(dr);
                    add_histogram_avg(s,i,dyn_ranges_time[flag][2*t],dyn_ranges_time[flag][2*t+1],dr);
                }
            }
        }
        eval_isoconf(t, flag);

        // the main evaluation for the isoconfigurational ensemble (log of mean displacement)
        //std::cout << "EVAL LOG(MD) " << t << std::endl; 
        reset_dyn(t);
        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) {
                for (int j=0; j<NI;j++) {
                    double dr = 0, dx;
                    for (int d=0; d<dim;d++) {
                        dx = xyz_data[i+s*N][d+dim*NT*j] - xyz_data[i+s*N][d+t*dim+dim*NT*j];
                        apply_pbc(dx);
                        dr += dx*dx;
                    }
                    dr = sqrt(dr);
                    //std::cout << dr << std::endl;
                    if (dr < 10e-6) dr = 10e-6;
                    dr = log10(dr);
                    //std::cout << dr << std::endl;
                    add_histogram_avg(s,i,dyn_ranges_time[flag+1][2*t],dyn_ranges_time[flag+1][2*t+1],dr);
                }
            }
        }
        //std::cout << "FINISHED LOG(MD) " << t << std::endl; 
        eval_isoconf(t, flag+1);

    }

    // write results
    print_isoconf(flag);
    print_isoconf(flag+1);

}