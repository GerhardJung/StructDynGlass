#ifndef STRUCT_BASE_H
#define STRUCT_BASE_H

#define RC2 6.25
#define C0 0.04049023795
#define C2 -0.00970155098
#define C4 0.00062012616

void eval_struct_base();

void rescale_print_gr();
void calc_psi(int ** neighbors);
void eval_den_cg();

// help function
double determine_sigma(int iType, int jType);
double determine_epsilon(int iType, int jType);
double calc_epot(int iType, int jType,  double dist2);

#endif