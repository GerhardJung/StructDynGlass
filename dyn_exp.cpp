
#include "dyn_exp.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_exp(){

    // loop over time
    for (int t=1; t<NT; t++) {
        std::cout << "EVAL EXP " << t << std::endl; 
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
                    double C_loc= exp(-dr*dr*exp_scale4i);
                    add_histogram_avg(s,i,dyn_ranges[exp_flag][0],dyn_ranges[exp_flag][1],C_loc);
                }
            }
        }

        // the main evaluation for the isoconfigurational ensemble
        eval_isoconf(t, exp_flag);

    }

    // write results
    print_isoconf(exp_flag,"EXP");


}

