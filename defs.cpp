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
int NDynTotal=0;               // number of final dynamical descriptors
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
double dyn_rearrange_threshold; //threshold to differentiate between active and inactive
int dyn_rearrange_mode=0;
double * save_pat;
double * save_pat_traj;

int NStruct;          // number of strctural observables to be analyzed
int NStructTotal=0;
std::string * StructNames; 

int struct_base_flag=-1;       // flag for the very basic structural descriptors
int NHistoGr;
double rcut2;            // cutoff for neighbor search
int lmax=10;
int lmin=1;

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
std::string struct_filion_file;
double ** struct_filion_descriptor_list;

int struct_ml_flag=-1;

int struct_gnn_flag=-1;

int struct_voronoi_flag=-1;

int struct_read_flag=-1;
int struct_read_Ntotal=0;
std::string struct_read_file;

// DATA
int N;                      // number of particles
double boxL;
int dim;
std::string model;           // name of the potential model (for epot and hessian)
int NTYPE;                   // number of different types
double * type_cutoff;        // cutoffs to deal with polydisperse samples
int NT;
int * type_data;
double * dia_data;             // particle diameters
double ** xyz_data;
double ** xyz_inherent_data;
int * time_data;
int * NPerType;

double * global_properties;

// DYN
double ** dyn_hist_data;
double *dyn_avg;
double * dyn_sys_avg;
double *dyn_avg2;
double **dyn_pred;
double ** dyn_avg_save;
double ** dyn_avg_save_cg;

//STRUCT
double **struct_local;

double **struct_base_gr;

double **struct_local_filion;
double **struct_filion_classifiers_thermal;
double **struct_filion_classifiers_inherent;
double ** struct_mean_var; 

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
    dia_data = dvector(0,N*NS-1);
    NPerType = ivector(0,NTYPE-1);
    xyz_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    xyz_inherent_data = dmatrix(0,N*NS-1,0,NI*NT*dim-1);
    time_data = ivector(0,NT-1);


    global_properties = dvector(0,NTYPE*5-1);

    // allocate dyn data
    dyn_hist_data = dmatrix(0,N*NS-1,0,NHisto-1);
    dyn_avg = dvector(0,N*NS-1);
    dyn_sys_avg = dvector(0,NS*NI*NTYPE-1);
    dyn_avg2 = dvector(0,N*NS-1);
    if (NDynTotal>0) dyn_pred = dmatrix(0,NT*NDynTotal-1,0,6*NTYPE-1);
    if (NDynTotal>0) dyn_avg_save = dmatrix(0,N*NS-1,0,NDynTotal*(NT+1)-1);
    if (NDynTotal>0) dyn_avg_save_cg = dmatrix(0,N*NS-1,0,NDynTotal*(NT+1)*5-1);

    // allocate struct data
    if (NStructTotal>0) struct_local = dmatrix(0,NStructTotal*NCG-1,0,N*NS-1);
    if (struct_filion_flag>=0) struct_local_filion = dmatrix(0,19*NCG-1,0,N*NS-1);

    struct_base_gr = dmatrix(0,NTYPE*NTYPE,0,NHistoGr-1);
    if(struct_soft_modes_flag>=0 || rp_flag>=0) hessian = dmatrix(0,N*N-1,0,dim*dim-1);
    if(struct_soft_modes_flag>=0) {
        hessian_evectors = dmatrix(0,N*dim-1,0,N*dim-1);
        hessian_evalues = dvector(0,NS*N*dim-1);
        sm_histograms = dmatrix(0,NHistoSM-1,0,2);
        participation_ratio = dvector(0,NS*N*dim-1);
    }
    if (struct_filion_flag >= 0) {
        struct_filion_descriptor_list = dmatrix(0,500,0,2);
    }

    // allocate dyn-struct correlation data and histogramms
    if (NDynTotal>0) dyn_hist_iso = dmatrix(0,NT*NDynTotal-1,0,NTYPE*NHisto-1);
    if (NDynTotal>0) dyn_hist_val = dmatrix(0,NT*NDynTotal-1,0,NTYPE*NHisto-1);
    if (NDynTotal>0) dyn_ranges_time = dmatrix(0,NDynTotal-1,0,2*NT-1);
    if (NStructTotal>0) struct_ranges = dmatrix(0,NCG*NStructTotal-1,0,1);
    if (NStructTotal>0) struct_hist = dmatrix(0,NCG*NStructTotal-1,0,NTYPE*NHistoStruct-1);
    dyn_struct_hist_iso = dmatrix(0,NTYPE-1,0,NHisto*NHistoStruct-1);
    dyn_struct_hist_val = dmatrix(0,NTYPE-1,0,NHisto*NHistoStruct-1);
    if (NDynTotal>0) if (NStructTotal>0) dyn_struct_pred = dmatrix(0,NCG*NStructTotal*NT*NDynTotal-1,0,6*NTYPE-1);


    // initialize data
    for (int i = 0; i < NS*N; i++) {
        type_data[i] = 0;
        for (int j = 0; j < NI*NT*dim; j++) {
            xyz_data[i][j] = 0.0;
            xyz_inherent_data[i][j] = 0.0;
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
        for (int k=0; k<5; k++) {
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

    for (int t=0; t<NT*NDynTotal; t++) {
        for (int j=0; j<6*NTYPE; j++) dyn_pred[t][j] = 0.0;
        for (int j=0; j<NTYPE*NHisto; j++) {
            dyn_hist_iso[t][j] = 0.0;
            dyn_hist_val[t][j] = 0.0;
        }
    }
    for (int j = 0; j < NCG*NStructTotal*NDynTotal*NT; j++) {
        for (int k = 0; k < 6*NTYPE; k++) dyn_struct_pred[j][k] = 0.0;
    }

}