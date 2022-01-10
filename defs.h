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

#include "nrutil.h"
#include "read_write_lammps.h"

#define NTYPE 3
#define NCG 4

#define EPS 0.000000001

// OPTIONS
extern std::string lammpsIn;       // lammps xyz output folder for readin
extern std::string folderOut; // folder to output isoconfigurational data
extern int CnfStart, CnfStep;      // to open output files
extern double timestep;            // size of the timestep
extern int NS;                     // number of different inital structures
extern int NI;                     // number of different isoconfigurational trajectories
extern int NHisto;                 // number of bin in the histograms (dynamical observables)
extern int NHistoStruct;           // number of bin in the histograms (structural observables)

extern int NDyn;                   // number of dynamical observables to be analyzed
extern std::string * DynNames; 

extern int bb_flag;                // flags for dynamical observables (bond-breaking)
extern double rcuti2;
extern double rcuto2;              //inner and outer cutoff for bond breaking

extern int exp_flag;               // flags for dynamical observables (exponential)
extern double exp_scale4i;         // length scale for exponential decay

extern int isf_flag;               // flags for dynamical observables (isf)
extern double qisf;                // length scale for isf

extern int msd_flag;            // flag for dynamical variables (msd)

extern int NStruct;          // number of strctural observables to be analyzed
extern int NStructTotal;
extern std::string * StructNames; 

extern int struct_base_flag;       // flag for the very basic structural descriptors
extern int NHistoGr;
extern double rcut2;            // cutoff for neighbor search

extern int struct_soft_modes_flag;       // flag for structural descriptors connected to soft modes
extern double ** hessian;           // hessian matrix for soft mode analysis

// DATA
extern int N;                      // number of particles
extern double boxL;                // box size
extern int dim;                    // number of dimensions
extern int NT;                     // number of timesteps
extern int * type_data;
extern double ** xyz_data;
extern double ** xyz_inherent_data;
extern int * time_data;
extern int * NPerType;

extern double * global_properties;

//DYN
extern double **dyn_hist_data;
extern double *dyn_avg;
extern double *dyn_avg2;
extern double **dyn_pred;
extern double **dyn_avg_save;


//STRUCT
extern double **struct_local;

extern double **struct_base_gr;
extern double **struct_base_rad_classifier;
extern double **struct_base_ang_classifier;

// DYN STRUCT CORRELATION, HISTOGRAMMS
extern double **dyn_ranges;
extern double **dyn_hist_iso;
extern double **dyn_hist_val;
extern double **struct_ranges;
extern double **struct_hist;
extern double **dyn_struct_hist_iso;
extern double **dyn_struct_hist_val;
extern double **dyn_struct_pred;

// Functions
void allocate_storage();



#endif