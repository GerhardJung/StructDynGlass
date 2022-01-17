/* Implement structural descriptors for the machine-learning technique as described in arXiv:2105.05921v1 */

#include "struct_filion.h"
#include "defs.h"
#include "pbc.h"

void eval_struct_filion(){

    std::cout << "EVAL STRUCT FILION " << std::endl; 

    int iType,jType;
    for (int s=0; s<NS;s++) { // loop over structures

        for (int i=0; i<N;i++) { // loop over particles
            //iType = type_data[i+s*N];

            for (int j=0; j<N;j++) { // loop over particle pairs
                jType = type_data[j+s*N];
                double dx[dim];
                for (int d=0; d<dim;d++) {
                    dx[d] = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx[d]);
                }

                // eval radial descriptors
                eval_radial(i+s*N,jType,dx);

                // eval angular descriptors
                eval_angular(i+s*N,dx);
                
            } 

        }
    }

    if (struct_filion_mode==0) {
        write_descriptors_csv();
    } else {
        read_eval_struct_filion_ml();
    }

}


void eval_radial(int i,int jType,double * dx){



}

void eval_angular(int i,double * dx){

}
                

void write_descriptors_csv(){

}


void read_eval_struct_filion_ml(){

}