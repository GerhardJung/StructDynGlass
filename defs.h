#ifndef DEFS_H
#define DEFS_H

// libraries
#include <iostream>
#include <math.h>       /* sqrt */
 
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QDebug>
#include <QString>
#include <chrono>

#include "nrutil.h"
#include "read_write_lammps.h"

#define EPS 0.000000001
#define NDim 2
//#define USE_RELATIVE  
#define EQUIBB     

// OPTIONS
extern std::string lammpsIn;       // lammps xyz output folder for readin
extern std::string folderOut; // folder to output isoconfigurational data
extern long CnfStart, CnfStep;      // to open output files
extern double timestep;            // size of the timestep
extern int NS;                     // number of different inital structures
extern int NI;                     // number of different isoconfigurational trajectories
extern int NHisto;                 // number of bin in the histograms (dynamical observables)
extern int NHistoStruct;           // number of bin in the histograms (structural observables)
extern int inherent;
extern int noinherent;
extern int NCG;
extern double * RCG;

extern int NDyn;                   // number of dynamical observables to be analyzed
extern int NDynTotal;               // number of final dynamical descriptors
extern std::string * DynNames; 

extern int bb_flag;                // flags for dynamical observables (bond-breaking)
extern double rcuti2;
extern double rcuto2;              //inner and outer cutoff for bond breaking

extern int exp_flag;               // flags for dynamical observables (exponential)
extern double exp_scale4i;         // length scale for exponential decay

extern int isf_flag;               // flags for dynamical observables (isf)
extern double qisf;                // length scale for isf

extern int msd_flag;            // flag for dynamical variables (msd)
extern double overlap_cut;

extern int rp_flag;             // flag for dynamical variables (strctural rearrangements as described by Patinet)
extern double dyn_rearrange_threshold; //threshold to differentiate between active and inactive
extern int dyn_rearrange_mode;
extern double * save_pat;
extern double * save_pat_traj;

extern int NStruct;          // number of strctural observables to be analyzed
extern int NStructTotal;
extern std::string * StructNames; 

extern int struct_base_flag;       // flag for the very basic structural descriptors
extern int NHistoGr;
extern double rcut2;            // cutoff for neighbor search
extern int lmax;
extern int lmin;

extern int struct_soft_modes_flag;       // flag for structural descriptors connected to soft modes
extern int NHistoSM;                    // number of histogram bins used for density of states evaluation
extern double Pcut;                 // cutoff for quasilocalized modes
extern double ** hessian;           // hessian matrix for soft mode analysis
extern double ** hessian_evectors;           // eigenvectors of the hessian matrix
extern double * hessian_evalues;      // eigenvalues of the hessian matrix
extern double ** sm_histograms;
extern double * participation_ratio;
extern int modeSM;

extern int struct_filion_flag;
extern std::string struct_filion_file;
extern double ** struct_filion_descriptor_list;

extern int struct_ml_flag;
extern int Ndes;

extern int struct_gnn_flag;

extern int struct_voronoi_flag;

extern int struct_read_flag;
extern int struct_read_Ntotal;
extern std::string struct_read_file;

// DATA
extern int N;                      // number of particles
extern double boxL;                // box size
extern std::string model;           // name of the potential model (for epot and hessian)
extern int NTYPE;                   // number of different types
extern double * type_cutoff;        // cutoffs to deal with polydisperse samples
extern int NT;                     // number of timesteps
extern int * type_data;             // particle types
extern double * dia_data;             // particle diameters
extern double ** xyz_data;
extern double ** xyz_inherent_data;
extern int * time_data;
extern int * NPerType;

extern double * global_properties;

//DYN
extern double **dyn_hist_data;
extern double *dyn_avg;
extern double *dyn_sys_avg;
extern double *dyn_avg2;
extern double **dyn_pred;
extern double **dyn_avg_save;
extern double **dyn_avg_save_cg;
extern int Nequi;
extern double * equiBB;

//STRUCT
extern double **struct_local;

extern double **struct_base_gr;

extern double **struct_filion_classifiers_thermal;
extern double **struct_filion_classifiers_inherent;
extern double **struct_local_ml;
extern double ** struct_mean_var; 

// DYN STRUCT CORRELATION, HISTOGRAMMS
extern double **dyn_ranges;
extern double **dyn_ranges_time;
extern double **dyn_hist_iso;
extern double **dyn_hist_val;
extern double **struct_ranges;
extern double **struct_hist;
extern double **dyn_struct_hist_iso;
extern double **dyn_struct_hist_val;
extern double **dyn_struct_pred;

extern double hist_lower_time;
extern double hist_upper_time;

// Functions
void allocate_storage();



#endif