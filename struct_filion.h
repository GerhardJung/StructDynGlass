#ifndef STRUCT_FILION_H
#define STRUCT_FILION_H

#define LMAX 10
#define CG_NMAX 2
#define CG_RC 2.3

void eval_struct_filion();

void init_descriptors();

void eval_radial(int i,int jType,double * dx, double dr);
void eval_angular(int i,double * dx, double dr);                

void normalize_cg_descriptors();

void write_descriptors_csv();

#endif