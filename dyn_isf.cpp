
#include "dyn_isf.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_isf(){

    for (int d=0; d<dim;d++) {
        // loop over time
        for (int t=1; t<NT; t++) {
            std::cout << "EVAL ISF DIM " << d << " " << t << std::endl; 
            reset_dyn(t);

            for (int s=0; s<NS;s++) { // loop over structures
                for (int i=0; i<N;i++) {
                    for (int j=0; j<NI;j++) {
                        double dx = xyz_data[i+s*N][d+dim*NT*j] - xyz_data[i+s*N][d+t*dim+dim*NT*j];
                        apply_pbc(dx);
                        double C_loc= cos(qisf*dx);
                        add_histogram_avg(s,i,dyn_ranges[isf_flag][0],dyn_ranges[isf_flag][1],C_loc);
                    }
                }
            }

            // the main evaluation for the isoconfigurational ensemble
            eval_isoconf(t,isf_flag);
        }
    }

    // write results
    print_isoconf(isf_flag,"ISF");
}

