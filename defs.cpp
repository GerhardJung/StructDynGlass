#include "defs.h"
#include "nrutil.h"

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

// DATA
int N;                      // number of particles
double boxL;
int dim;
int NT;
int * type_data;
double ** xyz_data;
int * time_data;

// DYN
double ** dyn_data;
double ** dyn_bb_data;

void allocate_storage(){

    // allocate type and xyz data
    type_data = ivector(0,N*NS-1);
    xyz_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    time_data = ivector(0,NT-1);

    // allocate dyn data
    dyn_data = dmatrix(0,N*NS-1,0,NHisto-1);
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
            dyn_data[i][j] = 0.0;
        }
    }

}