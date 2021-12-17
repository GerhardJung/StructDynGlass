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

// OPTIONS
extern std::string lammpsIn;       // lammps xyz output folder for readin
extern std::string folderOut; // folder to output isoconfigurational data
extern int CnfStart, CnfStep;      // to open output files
extern double timestep;            // size of the timestep
extern int NS;                     // number of different inital structures
extern int NI;                     // number of different isoconfigurational trajectories
extern int NHisto;                 // number of histograms
extern int NDyn;                   // number of dynamical observables to be analyzed

extern int bb_flag;                // flags for dynamical observables (bond-breaking)
extern double rcuti2;
extern double rcuto2;              //inner and outer cutoff for bond breaking
extern double bb_hist_lower;
extern double bb_hist_upper;

extern int exp_flag;               // flags for dynamical observables (exponential)
extern double exp_scale4i;         // length scale for exponential decay
extern double exp_hist_lower;
extern double exp_hist_upper;

extern int isf_flag;               // flags for dynamical observables (isf)
extern double qisf;                // length scale for isf
extern double isf_hist_lower;
extern double isf_hist_upper;

extern int NStruct;          // number of strctural observables to be analyzed

extern int struct_base_flag;       // flag for the very basic structural descriptors
extern int NHistoGr;

// DATA
extern int N;                      // number of particles
extern double boxL;                // box size
extern int dim;                    // number of dimensions
extern int NT;                     // number of timesteps
extern int * type_data;
extern double ** xyz_data;
extern int * time_data;
extern int * NPerType;

//DYN
extern double **dyn_hist_data;
extern double *dyn_avg;
extern double *dyn_avg2;
extern double **dyn_pred;
extern double **dyn_hist_iso;
extern double **dyn_hist_val;

extern double **dyn_bb_avg;
extern double **dyn_exp_avg;
extern double **dyn_isf_avg;

//STRCUT
extern double **struct_base_gr;
extern double **struct_base_theta5;
extern double **struct_base_theta6;
extern double **struct_base_rad_classifier;
extern double **struct_base_ang_classifier;
extern double **struct_base_local_den;
extern double **struct_base_local_epot;
extern double **struct_base_local_theta;
extern double **struct_base_local_psi5;
extern double **struct_base_local_psi6;

// Functions
void allocate_storage();

#endif