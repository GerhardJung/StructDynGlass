#ifndef EVAL_ISOCONF_H
#define EVAL_ISOCONF_H

#include "defs.h"

void reset_dyn_hist();

void add_histogram(int s, int i, double val);

void norm_histogram();

void eval_isoconf(int t);

#endif