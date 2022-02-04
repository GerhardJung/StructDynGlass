#ifndef STRUCT_BASE_H
#define STRUCT_BASE_H

#define RC2LJ 6.25
#define C0LJ 0.04049023795
#define C2LJ -0.00970155098
#define C4LJ 0.00062012616

#define RC2POLY 1.5625
#define C0POLY -1.92414534861
#define C2POLY 2.11106232533
#define C4POLY -0.5910974510

void eval_struct_base();

void rescale_print_gr();
void calc_psi(int ** neighbors);
void eval_den_cg();

// help function
double determine_sigma(int i, int j);
double determine_epsilon(int iType, int jType);
double calc_epot(int i, int j,  double dist2);

#endif