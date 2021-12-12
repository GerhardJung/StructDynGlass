
#include "pbc.h"
#include "defs.h"


void apply_pbc_global(){
    // apply pbc
    for (int i = 0; i < NS*N; i++) {
        type_data[i] = 0;
        for (int j = 0; j < NI*NT*dim; j++) {
            apply_pbc(xyz_data[i][j]);
        }
    }

}

void apply_pbc(double &x){
    while (x < boxL/2.0) x += boxL; 
    while (x > boxL/2.0) x -= boxL; 
}