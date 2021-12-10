
// OPTIONS
std::string lammpsIn;       // lammps xyz output file for readin
std::string xyzOut;         // extended xyz file to output isoconfigurational data
int CnfStart, CnfStep;      // to open output files
int NS;                     // number of different inital structures
int NI;                     // number of different isoconfigurational trajectories
int NDyn;                   // number of dynamical observables to be analyzed
int bb_flag;                // flags for dynamical observables