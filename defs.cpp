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
std::string * DynNames; 

int bb_flag=-1;                // flags for dynamical observables
double rcuti2;
double rcuto2;              //inner and outer cutoff for bond breaking

int exp_flag=-1;               // flags for dynamical observables (exponential)
double exp_scale4i;         // length scale for exponential decay

int isf_flag=-1;               // flags for dynamical observables (isf)
double qisf;                // length scale for isf

int msd_flag=-1;            // flag for dynamical variables (msd)

int rp_flag=-1;             // flag for dynamical variables (strctural rearrangements as described by Patinet)

int NStruct;          // number of strctural observables to be analyzed
int NStructTotal=0;
std::string * StructNames; 

int struct_base_flag=-1;       // flag for the very basic structural descriptors
int NHistoGr;
double rcut2;            // cutoff for neighbor search

int struct_soft_modes_flag=-1;       // flag for structural descriptors connected to soft modes
int NHistoSM;                    // number of histogram bins used for density of states evaluation
double Pcut;                 // cutoff for quasilocalized modes
double ** hessian;           // hessian matrix for soft mode analysis
double ** hessian_evectors;           // eigenvectors of the hessian matrix
double * hessian_evalues;      // eigenvalues of the hessian matrix
double ** sm_histograms;
int modeSM;
double * participation_ratio;

int struct_filion_flag=-1;
int struct_filion_mode;

// DATA
int N;                      // number of particles
double boxL;
int dim;
int NT;
int * type_data;
double ** xyz_data;
double ** xyz_inherent_data;
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
double **struct_local;

double **struct_base_gr;
double **struct_base_rad_classifier;
double **struct_base_ang_classifier;

// DYN STRUCT CORRELATION, HISTOGRAMMS
double **dyn_ranges;
double **dyn_ranges_time;
double **dyn_hist_iso;
double **dyn_hist_val;
double **struct_ranges;
double **struct_hist;
double **dyn_struct_hist_iso;
double **dyn_struct_hist_val;
double **dyn_struct_pred;

void allocate_storage(){

    // allocate type and xyz data
    type_data = ivector(0,N*NS-1);
    xyz_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    xyz_inherent_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    time_data = ivector(0,NT-1);

    global_properties = dvector(0,NTYPE*2-1);

    // allocate dyn data
    dyn_hist_data = dmatrix(0,N*NS-1,0,NHisto-1);
    dyn_avg = dvector(0,N*NS-1);
    dyn_avg2 = dvector(0,N*NS-1);
    dyn_pred = dmatrix(0,NT-1,0,6*NTYPE-1);
    dyn_avg_save = dmatrix(0,N*NS-1,0,NDyn*NT-1);

    // allocate struct data
    struct_local = dmatrix(0,NStructTotal*NCG-1,0,N*NS-1);

    struct_base_gr = dmatrix(0,NTYPE*NTYPE,0,NHistoGr-1);
    hessian = dmatrix(0,NS*N*N-1,0,dim*dim-1);
    hessian_evectors = dmatrix(0,NS*N*dim-1,0,N*dim-1);
    hessian_evalues = dvector(0,NS*N*dim-1);
    sm_histograms = dmatrix(0,NHistoSM-1,0,2);
    participation_ratio = dvector(0,NS*N*dim-1);

    // allocate dyn-struct correlation data and histogramms
    dyn_hist_iso = dmatrix(0,NT-1,0,NTYPE*NHisto-1);
    dyn_hist_val = dmatrix(0,NT-1,0,NTYPE*NHisto-1);
    dyn_ranges_time = dmatrix(0,NDyn-1,0,NT*2-1);
    struct_ranges = dmatrix(0,NCG*NStructTotal-1,0,1);
    struct_hist = dmatrix(0,NCG*NStructTotal-1,0,NTYPE*NHistoStruct-1);
    dyn_struct_hist_iso = dmatrix(0,NTYPE-1,0,NHisto*NHistoStruct-1);
    dyn_struct_hist_val = dmatrix(0,NTYPE-1,0,NHisto*NHistoStruct-1);
    dyn_struct_pred = dmatrix(0,NCG*NStructTotal*NT-1,0,6*NTYPE-1);


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

        for (int j = 0; j < NCG*NStructTotal; j++) {
            struct_local[j][i] = 0.0;
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