#ifndef GLOBAL_H
#define GLOBAL_H

void calc_print_global();

void calc_bonds_histograms_structure();

// help function that writes structural information with index i into array 
void struct_array(int index, int c, double * array);
void copy_array(double ** array_in, int index, double * array_out);

#endif