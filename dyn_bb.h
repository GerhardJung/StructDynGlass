#ifndef DYN_BB_H
#define DYN_BB_H

#define N_NEIGH_MAX 25

void eval_bb();

// Help functions
void findneighbors(int rcut2, int ** neighbors);
void checkneighbors(int s, int i, int j, int t, int &n0, int &nt, int ** neighbors);

#endif