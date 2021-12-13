#include "defs.h"

//OPTIONS
std::string lammpsIn;       // lammps xyz output file for readin
std::string xyzOut;         // extended xyz file to output isoconfigurational data
int CnfStart, CnfStep;      // to open output files
double timestep;
int NS;                     // number of different inital structures
int NI;                     // number of different isoconfigurational trajectories
int NHisto;
int NDyn;                   // number of dynamical observables to be analyzed

int bb_flag;                // flags for dynamical observables
double rcuti2;
double rcuto2;              //inner and outer cutoff for bond breaking

// DATA
int N;                      // number of particles
double boxL;
int dim;
int NT;
int * type_data;
double ** xyz_data;
int * time_data;
int * NPerType = new int[NTYPE];

// DYN
double ** dyn_hist_data;
double **dyn_pred;
double **dyn_hist_iso;
double **dyn_hist_val;

double ** dyn_bb_data;

void allocate_storage(){

    // allocate type and xyz data
    type_data = ivector(0,N*NS-1);
    xyz_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    time_data = ivector(0,NT-1);

    // allocate dyn data
    dyn_hist_data = dmatrix(0,N*NS-1,0,NHisto-1);
    dyn_pred = dmatrix(0,NT-1,0,5*NTYPE-1);
    dyn_hist_iso = dmatrix(0,NT-1,0,NTYPE*NHisto-1);
    dyn_hist_val = dmatrix(0,NT-1,0,NTYPE*NHisto-1);

    dyn_bb_data = dmatrix(0,N*NS-1,0,NT-1);

    // initialize data
    for (int i = 0; i < NS*N; i++) {
        type_data[i] = 0;
        for (int j = 0; j < NI*NT*dim; j++) {
            xyz_data[i][j] = 0.0;
        }
    }

    for (int i = 0; i < NS*N; i++) {
        for (int j = 0; j < NT; j++) {
            dyn_bb_data[i][j] = 0.0;
        }
        for (int j = 0; j < NHisto; j++) {
            dyn_hist_data[i][j] = 0.0;
        }
    }

    for (int t = 0; t < NT; t++) {
        for (int j = 0; j < NTYPE*NHisto; j++) {
            dyn_hist_iso[t][j] = 0.0;
            dyn_hist_val[t][j] = 0.0;
        }
        for (int j = 0; j < 5*NTYPE; j++) {
            dyn_pred[t][j] = 0.0;
        }
    }

    for (int j = 0; j < NTYPE; j++) {
        NPerType[j] = 0;
    }

}