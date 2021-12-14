
#include "dyn_isf.h"
#include "defs.h"
#include "pbc.h"
#include "eval_isoconf.h"

void eval_isf(){

    for (int d=0; d<dim;d++) {
        // loop over time
        for (int t=1; t<NT; t++) {
            std::cout << "EVAL ISF DIM " << d << " " << t << std::endl; 
            reset_dyn();

            for (int s=0; s<NS;s++) { // loop over structures
                for (int i=0; i<N;i++) {
                    for (int j=0; j<NI;j++) {
                        double dx = xyz_data[i+s*N][d+dim*NT*j] - xyz_data[i+s*N][d+t*dim+dim*NT*j];
                        apply_pbc(dx);
                        double C_loc= cos(qisf*dx);
                        add_histogram_avg(s,i,isf_hist_lower,isf_hist_upper,C_loc);
                    }
                }
            }

            // the main evaluation for the isoconfigurational ensemble
            if(d==0) eval_isoconf(t, "ISFX");
            if(d==1) eval_isoconf(t, "ISFY");
            if(d==2) eval_isoconf(t, "ISFZ");

            // write results
            if(d==0) print_isoconf("ISFX");
            if(d==1) print_isoconf("ISFY");
            if(d==2) print_isoconf("ISFZ");

        }
    }
}

