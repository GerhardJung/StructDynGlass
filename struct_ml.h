#ifndef STRUCT_ML_H
#define STRUCT_ML_H

// cell lists
extern int ** cell_list_index;
extern int **cell_list;
extern double rc;			// constants for cell list calculation
extern int Ncell;
extern int Nmax;

void eval_struct_ml();
void eval_den_cg_ml();
void write_descriptors_csv();
void calc_voronoi(int s);

void write_descriptors_csv_dyn();

void create_cell_lists();

#endif