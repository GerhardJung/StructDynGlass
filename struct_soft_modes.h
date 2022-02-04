#ifndef STRUCT_SOFT_MODES_H
#define STRUCT_SOFT_MODES_H

#define NLOW 30

void eval_struct_soft_modes();

// help function
void calc_2Depot(int i, int j, int t, int iso, double * result);
void calc_3Depot(int i, int j, int k, double * result);
double calc_epot_tot(int s);

#endif