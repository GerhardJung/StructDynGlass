#ifndef EVAL_ISOCONF_H
#define EVAL_ISOCONF_H

#include "defs.h"

void reset_dyn(int t);
void add_histogram_avg(int s, int i, double hist_lower, double hist_upper, double val);
void norm_histogram();

void eval_isoconf(int t, int loc);

void calc_histograms_dynamics(int t, int loc);
void eval_information_theory_dynamics(int t);

void calc_histograms_information_theory(double * struct_array, int loc_dyn,int loc_struct, int c);
void eval_information_theory_correlation(int t,int loc_struct,int c);
void eval_pearson_spearman_correlation(int t,double * struct_array,int loc_struct,int c);

void print_isoconf(int loc);

#endif