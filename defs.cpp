#include "defs.h"

//OPTIONS
std::string lammpsIn;       // lammps xyz output folder for readin
std::string folderOut;         // folder to output isoconfigurational data
int CnfStart, CnfStep;      // values to open output files
double timestep;
int NS;                     // number of different inital structures
int NI;                     // number of different isoconfigurational trajectories
int NHisto;
int NHistoStruct;           // number of bin in the histograms (structural observables)

int NDyn;                   // number of dynamical observables to be analyzed

int bb_flag=-1;                // flags for dynamical observables
double rcuti2;
double rcuto2;              //inner and outer cutoff for bond breaking
double bb_hist_lower;
double bb_hist_upper;

int exp_flag=-1;               // flags for dynamical observables (exponential)
double exp_scale4i;         // length scale for exponential decay
double exp_hist_lower;
double exp_hist_upper;

int isf_flag=-1;               // flags for dynamical observables (isf)
double qisf;                // length scale for isf
double isf_hist_lower;
double isf_hist_upper;

int NStruct;          // number of strctural observables to be analyzed
int NStructTotal;

int struct_base_flag=-1;       // flag for the very basic structural descriptors
int NHistoGr;
double rcut2;            // cutoff for neighbor search
int bo_Nneigh;    // how many neighbors are considered for bo average

// DATA
int N;                      // number of particles
double boxL;
int dim;
int NT;
int * type_data;
double ** xyz_data;
int * time_data;
int * NPerType = new int[NTYPE];

double * global_properties;

// DYN
double ** dyn_hist_data;
double *dyn_avg;
double *dyn_avg2;
double **dyn_pred;
double ** dyn_avg_save;

//STRUCT
double **struct_base_gr;

double **struct_base_rad_classifier;
double **struct_base_ang_classifier;
double **struct_base_local_den;
double **struct_base_local_epot;
double **struct_base_local_theta_tanaka;
double **struct_base_local_psi;
double **struct_base_local_bo;

// DYN STRUCT CORRELATION, HISTOGRAMMS
int **dyn_ranges;
double **dyn_hist_iso;
double **dyn_hist_val;
int **struct_ranges;
double **struct_hist;
double **dyn_struct_hist_iso;
double **dyn_struct_hist_val;
double **dyn_struct_pred;

void allocate_storage(){

    // allocate type and xyz data
    type_data = ivector(0,N*NS-1);
    xyz_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    time_data = ivector(0,NT-1);

    global_properties = dvector(0,NTYPE*2-1);

    // allocate dyn data
    dyn_hist_data = dmatrix(0,N*NS-1,0,NHisto-1);
    dyn_avg = dvector(0,N*NS-1);
    dyn_avg2 = dvector(0,N*NS-1);
    dyn_pred = dmatrix(0,NT-1,0,5*NTYPE-1);

    dyn_avg_save = dmatrix(0,NDyn*N*NS-1,0,NT-1);

    // allocate struct data
    struct_base_gr = dmatrix(0,NTYPE*NTYPE,0,NHistoGr-1);

    struct_base_local_psi = dmatrix(0,N*NS-1,0,2*NCG-1);
    //struct_base_local_bo = dmatrix(0,N*NS-1,0,2-1);
    struct_base_local_theta_tanaka = dmatrix(0,N*NS-1,0,NCG-1);
    struct_base_local_epot = dmatrix(0,N*NS-1,0,NCG-1);
    struct_base_local_den = dmatrix(0,N*NS-1,0,NCG-1);

    // allocate dyn-struct correlation data and histogramms
    dyn_hist_iso = dmatrix(0,NT-1,0,NTYPE*NHisto-1);
    dyn_hist_val = dmatrix(0,NT-1,0,NTYPE*NHisto-1);
    struct_ranges = imatrix(0,NCG*NStructTotal-1,0,1);
    struct_hist = dmatrix(0,NCG*NStructTotal-1,0,NTYPE*NHistoStruct-1);
    dyn_struct_hist_iso = dmatrix(0,NTYPE-1,0,NHisto*NHistoStruct-1);
    dyn_struct_hist_val = dmatrix(0,NTYPE-1,0,NHisto*NHistoStruct-1);
    dyn_struct_pred = dmatrix(0,NCG*NStructTotal*NT-1,0,4*NTYPE-1);


    // initialize data
    for (int i = 0; i < NS*N; i++) {
        type_data[i] = 0;
        for (int j = 0; j < NI*NT*dim; j++) {
            xyz_data[i][j] = 0.0;
        }
    }

    for (int i = 0; i < NS*N; i++) {
        dyn_avg2[i] = 0.0;
        dyn_avg[i] = 0.0;
        for (int j = 0; j < NHisto; j++) {
            dyn_hist_data[i][j] = 0.0;
        }

        for (int j = 0; j < NCG; j++) {
            struct_base_local_psi[i][2*j] = 0.0;
            struct_base_local_psi[i][2*j+1] = 0.0;
            //struct_base_local_bo[i][0] = 0.0;
            struct_base_local_theta_tanaka[i][j] = 0.0;
            struct_base_local_epot[i][j] = 0.0;
            struct_base_local_den[i][j] = 0.0;

        }
    }

    for (int j = 0; j < NTYPE; j++) {
        NPerType[j] = 0;
        for (int k=0; k<2; k++) {
            global_properties[k*NTYPE+j] = 0.0;
        }
    }

    for (int i = 0; i < NTYPE; i++) {
        for (int j = 0; j <= NTYPE; j++) {
            for (int k=0; k<NHistoGr; k++) {
                struct_base_gr[i*NTYPE+j][k]=0.0;
            }
        }
    }

}