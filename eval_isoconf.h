#ifndef EVAL_ISOCONF_H
#define EVAL_ISOCONF_H

#include "defs.h"

void reset_dyn();

void add_histogram_avg(int s, int i, double hist_lower, double hist_upper, double val);

void norm_histogram();

void eval_isoconf(int t, std::string dyn);

void print_isoconf(std::string dyn);

#endif