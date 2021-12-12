
#include "dyn_bb.h"
#include "defs.h"
#include "pbc.h"

#define N_NEIGH_MAX 15

int **neighbors;

void eval_bb(){
    neighbors = imatrix(0,NS*N-1,0,N_NEIGH_MAX-1);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            for (int d=0; d<N_NEIGH_MAX;d++) {
                neighbors[i+s*N][d] = -1;
            }
        }
    }
}



// Help functions
void findneighbors() {
    double dr, dx;
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            int ncount = 0;
            for (int j=0; j<N;j++) {
                dr = 0;
                for (int d=0; d<dim;d++) {
                    dx = xyz_data[i+s*N][0] - xyz_data[j+s*N][0];
                    apply_pbc(dx);
                    dr += dx*dx;
                }
                

                if (dr < rcuti2 && i != j ) {

                    neighbors[i+s*N][ncount] = j;
                    ncount ++;

                }
            }
        }
    }
}

void checkneighbors(int s, int i, int &n0, int &nt) {
    double dr, dx;
    n0 = 0;
    nt = 0;
    while (neighbors[i+s*N][n0] != -1 ) {
        dr = 0;
        for (int d=0; d<dim;d++) {
            dx = xyz_data[i+s*N][0] - xyz_data[neighbors[i+s*N][n0]+s*N][0];
            apply_pbc(dx);
            dr += dx*dx;
        }

      if (dr < rcuto2) {
         nt ++;
      }
      n0 ++;
   }
}