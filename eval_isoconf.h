#ifndef EVAL_ISOCONF_H
#define EVAL_ISOCONF_H

#include "defs.h"

void reset_dyn(int t);
void add_histogram_avg(int s, int i, int j,double* hist_cut, double val);
void norm_histogram();

void eval_isoconf(int t, int flag);

void calc_histograms_dynamics(int t, int flag);
void eval_information_theory_dynamics(int t);

void calc_histograms_information_theory(double * struct_array, int loc_dyn,int loc_struct, int c);
void eval_information_theory_correlation(int t,int loc_struct,int c);
void eval_pearson_spearman_correlation(int t,double * struct_array,int loc_struct,int flag, int c);

void print_isoconf(int flag);
void print_traj(double * save_dyn);

#endif