
#include "pbc.h"
#include "defs.h"


void apply_pbc_global(){
    // apply pbc
    for (int i = 0; i < NS*N; i++) {
        for (int j = 0; j < NI*NT*dim; j++) {
            apply_pbc(xyz_data[i][j]);
            apply_pbc(xyz_inherent_data[i][j]);
        }
    }

}

void apply_pbc(double &x){
    while (x < boxL/2.0) x += boxL; 
    while (x > boxL/2.0) x -= boxL; 
}