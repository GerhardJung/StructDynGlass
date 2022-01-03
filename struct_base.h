#ifndef STRUCT_BASE_H
#define STRUCT_BASE_H

void eval_struct_base();

void rescale_print_gr();
void calc_psi(int ** neighbors);
void eval_epot_den_cg();

// help function
double determine_sigma(int iType, int jType);
double determine_epsilon(int iType, int jType);
double calc_epot(int i, int j,  double dist);

#endif