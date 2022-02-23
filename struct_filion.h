#ifndef STRUCT_FILION_H
#define STRUCT_FILION_H

#define CG_NMAX 2
#define CG_RC 2.3

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

void eval_struct_filion();

void init_descriptors();

void eval_radial(int i,int jType,double * dx, double dr, double ** out);
void eval_angular(int i,double * dx, double dr, double * out);                

void collect_angular_descriptors(int i, double * struct_filion_save,  double ** struct_filion_classifiers);
void normalize_cg_descriptors(int struct_filion_mode,double ** struct_filion_classifiers);

void write_descriptors_csv_phys();
void write_descriptors_csv_dyn();

#endif