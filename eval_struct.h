#ifndef EVAL_STRUCT_H
#define EVAL_STRUCT_H

#include "defs.h"

void eval_struct(double * input_tensor,std::string input_name, int first);

void calc_bonds_histograms_structure();

// help function
void calc_bonds(double * input, double * output);

#endif