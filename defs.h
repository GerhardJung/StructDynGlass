#ifndef DEFS_H
#define DEFS_H

// libraries
#include <iostream>
 
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QDebug>
#include <QString>

// OPTIONS
extern std::string lammpsIn;       // lammps xyz output file for readin
extern std::string xyzOut;         // extended xyz file to output isoconfigurational data
extern int CnfStart, CnfStep;      // to open output files
extern double timestep;            // size of the timestep
extern int NS;                     // number of different inital structures
extern int NI;                     // number of different isoconfigurational trajectories
extern int NHisto;                 // number of histograms
extern int NDyn;                   // number of dynamical observables to be analyzed
extern int bb_flag;                // flags for dynamical observables

// DATA
extern int N;                      // number of particles
extern double boxL;                // box size
extern int dim;                    // number of dimensions
extern int NT;                     // number of timesteps
extern int * type_data;
extern double ** xyz_data;
extern int * time_data;

//DYN
extern double **dyn_data;
extern double **dyn_bb_data;

// Functions
void allocate_storage();

#endif