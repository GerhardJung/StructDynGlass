#include "struct_soft_modes.h"
#include "struct_base.h"
#include "defs.h"
#include "pbc.h"
#include <Eigen/Dense>

using namespace Eigen;

void eval_struct_soft_modes(){

    std::cout << "EVAL STRUCT SOFT MODES " << std::endl; 

    // calculate hessian
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles
            for (int j=0; j<N;j++) { // loop over particles
                calc_2Depot(i+s*N,j+s*N,hessian[s*N*N+i*N+j]);
            }
        }
    }

    // test harmonic
    /*int s=0;
    int i=1;
    int j=2;
    int d1=0;
    int d2=1;
    double shift = 0.00;
    double delta = 0.000001;
    xyz_inherent_data[i+s*N][d1] -= shift;
    xyz_inherent_data[i+s*N][d1] -= delta;
    double epot_im1 = calc_epot_tot(s);
    xyz_inherent_data[i+s*N][d1] += 2.0*delta;
    double epot_ip1 = calc_epot_tot(s);
    xyz_inherent_data[i+s*N][d1] -= delta;
    xyz_inherent_data[i+s*N][d1] += shift;

    //std::cout << epot_ip1 << " " << epot_im1 << " " <<  0.5/delta*(epot_ip1-epot_im1) << std::endl;

    // test hessian
    xyz_inherent_data[i+s*N][d1] -= delta;
    xyz_inherent_data[j+s*N][d2] -= delta;
    double epot_im1_jm1 = calc_epot_tot(s);
    xyz_inherent_data[i+s*N][d1] += 2.0*delta;
    double epot_ip1_jm1 = calc_epot_tot(s);
    xyz_inherent_data[j+s*N][d2] += 2.0*delta;
    double epot_ip1_jp1 = calc_epot_tot(s);
    xyz_inherent_data[i+s*N][d1] -= 2.0*delta;
    double epot_im1_jp1 = calc_epot_tot(s);
    xyz_inherent_data[i+s*N][d1] += delta;
    xyz_inherent_data[j+s*N][d2] -= delta;
    std::cout << hessian[s*N*N+i*N+j][d1*dim+d2] << " " << 0.25/delta/delta*(epot_ip1_jp1-epot_im1_jp1-epot_ip1_jm1+epot_im1_jm1) << std::endl;*/

    // calculate eigen decomposition
    MatrixXd Eigen_hessian = MatrixXd::Zero(N*dim,N*dim);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles
            for (int d1=0; d1<dim;d1++) {
                for (int j=0; j<N;j++) { // loop over particles
                    for (int d2=0; d2<dim;d2++) {
                        Eigen_hessian(i*d1,j*d2) = hessian[s*N*N+i*N+j][d1*dim+d2];
                    }
                }
            }
        }

        // the actual eigen decomposition
        EigenSolver<MatrixXd> Eigen_decomposition(Eigen_hessian);
        std::cout << "The largest eigenvalues of A is:" << std::endl << Eigen_decomposition.eigenvalues()[0] << std::endl;
        VectorXcd vT = Eigen_decomposition.eigenvectors().col(0);
        std::cout << "The corresponding eigenvector is:" << std::endl << vT << std::endl; 

        for (int k=0; k<N*dim;k++) { 
            hessian_evalues[s*N*dim+k] = Eigen_decomposition.eigenvalues()[N*dim-1-k].real();
            VectorXcd vT = Eigen_decomposition.eigenvectors().col(N*dim-1-k);
            for (int i=0; i<N*dim;i++) { 
                hessian_evectors[s*N*dim+k][i] = vT(i).real();
            }
        }
    }

    // calculate participation ratio for lowest frequency modes
    for (int s=0; s<NS;s++) {
        for (int i=0; i<N;i++) {
            struct_local[NCG*(struct_soft_modes_flag)][i+s*N] = 0.0;
            for (int k=0; k<NLOW;k++) {
                // generalize to dim dimensions!
                for (int di=0; di<dim;di++) struct_local[NCG*(struct_soft_modes_flag)][i+s*N] += hessian_evectors[s*N*dim+k][i*dim+di]*hessian_evectors[s*N*dim+k][i*dim+di];
            }   
        }
    }

    // calculate vibrality
    for (int s=0; s<NS;s++) {
        for (int i=0; i<N;i++) {
            struct_local[NCG*(struct_soft_modes_flag+1)][i+s*N] = 0.0;
            for (int k=0; k<N*dim;k++) {
                for (int di=0; di<dim;di++) struct_local[NCG*(struct_soft_modes_flag+1)][i+s*N] += hessian_evectors[s*N*dim+k][i*dim+di]*hessian_evectors[s*N*dim+k][i*dim+di];
                struct_local[NCG*(struct_soft_modes_flag+1)][i+s*N] /= hessian_evalues[s*N*dim+k]*hessian_evalues[s*N*dim+k];
            }   
        }
    }
}





// help functions

// calc hessian matrix
void calc_2Depot(int i, int j, double * result) {

    if (i==j) { // need to sum over all neighbors
        double result_loc[dim*dim];

        // calculate types
        int iType = type_data[i];
        int s = (i - i % N)/N;
        //std::cout << "s " << s << std::endl;
        for (int k=0; k<N;k++) { // loop over particle pairs
            int kType = type_data[k+s*N];
            double dr = 0.0, dx[dim];
            for (int d=0; d<dim;d++) {
                dx[d] = xyz_inherent_data[i][d] - xyz_inherent_data[k+s*N][d];
                apply_pbc(dx[d]);
                dr += dx[d]*dx[d];
            }

            if (dr < RC2 && i!=k ) {
                //std::cout << i << " " << k << std::endl;
                for (int d1=0; d1<dim;d1++) {
                    for (int d2=0; d2<dim;d2++) {
                        result_loc[d1*dim+d2] = 0.0;
                    }
                }
                double sigma = determine_sigma(iType, kType);
                double epsilon = determine_epsilon(iType, kType);
                double rij2i =  1.0/dr;
                double sigma2= sigma*sigma;
                double rij6i = sigma2*sigma2*sigma2*rij2i*rij2i*rij2i;
                double rij12i = rij6i*rij6i;
                double c2 = C2LJ / (sigma2);
                double c4 = C4LJ / (sigma2*sigma2);
                for (int d1=0; d1<dim;d1++) {
                    for (int d2=0; d2<dim;d2++) {
                        if (d1==d2) result_loc[d1*dim+d2] = 2.0*c2 + 8.0*c4*dx[d1]*dx[d1] + 4.0*c4*dr + 168.0*dx[d1]*dx[d1]*rij12i*rij2i*rij2i - 12.0*rij12i*rij2i - 48.0*dx[d1]*dx[d1]*rij6i*rij2i*rij2i + 6.0*rij6i*rij2i;
                        else {
                            if (d1 > d2) result_loc[d1*dim+d2] = result_loc[d2*dim+d1];
                            else result_loc[d1*dim+d2] = 8.0*c4*dx[d1]*dx[d2] + 168.0*dx[d1]*dx[d2]*rij12i*rij2i*rij2i - 48.0*dx[d1]*dx[d2]*rij6i*rij2i*rij2i;
                        }
                        result[d1*dim+d2] += 4.0*epsilon*result_loc[d1*dim+d2];
                    }
                }
            }
        }

    } else {
        // calculate types
        int iType = type_data[i];
        int jType = type_data[j];

        // calculate interaction parameters
        double sigma = determine_sigma(iType, jType);
        double epsilon = determine_epsilon(iType, jType);

        double dr = 0.0, dx[dim];
        for (int d=0; d<dim;d++) {
            dx[d] = xyz_inherent_data[i][d] - xyz_inherent_data[j][d];
            apply_pbc(dx[d]);
            dr += dx[d]*dx[d];
        }

        if (dr < RC2 ) {
            double result_loc[dim*dim];
            double rij2i =  1.0/dr;
            double sigma2= sigma*sigma;
            double rij6i = sigma2*sigma2*sigma2*rij2i*rij2i*rij2i;
            double rij12i = rij6i*rij6i;
            double c2 = C2LJ / (sigma2);
            double c4 = C4LJ / (sigma2*sigma2);
            for (int d1=0; d1<dim;d1++) {
                for (int d2=0; d2<dim;d2++) {
                    if (d1==d2) result_loc[d1*dim+d2] = 2.0*c2 + 8.0*c4*dx[d1]*dx[d1] + 4.0*c4*dr + 168.0*dx[d1]*dx[d1]*rij12i*rij2i*rij2i - 12.0*rij12i*rij2i - 48.0*dx[d1]*dx[d1]*rij6i*rij2i*rij2i + 6.0*rij6i*rij2i;
                    else {
                        if (d1 > d2) result_loc[d1*dim+d2] = result_loc[d2*dim+d1];
                        else result_loc[d1*dim+d2] = 8.0*c4*dx[d1]*dx[d2] + 168.0*dx[d1]*dx[d2]*rij12i*rij2i*rij2i - 48.0*dx[d1]*dx[d2]*rij6i*rij2i*rij2i;
                    }
                    result[d1*dim+d2] = - 4.0*epsilon*result_loc[d1*dim+d2];
                }
            }
        } 
    }
}

void calc_3Depot(int i, int j, int k, double * result) {
    /*double sigma = determine_sigma(iType, jType);
    double epsilon = determine_epsilon(iType, jType);

    if (dist2 < RC2 ) {
        double rij2 = dist2/(sigma*sigma);
        double rij4 = rij2*rij2;
        double rij6 = 1.0/(rij4*rij2);
        double rij12 = rij6*rij6;
         4.0*epsilon*(C0+C2*rij2+C4*rij4 - rij6 + rij12 ) ;
    } */
}

double calc_epot_tot(int s){
    double epot=0.0;
    for (int i=0; i<N;i++) { // loop over particles
        int iType = type_data[i+s*N];

        for (int j=0; j<N;j++) { // loop over particle pairs
            int jType = type_data[j+s*N];
            double dr = 0.0, dx;
            for (int d=0; d<dim;d++) {
                dx = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                apply_pbc(dx);
                dr += dx*dx;
            }

            // calculate epot
            if (i!=j) epot += 0.5*calc_epot(iType, jType, dr);
        }
    }
    return epot;
}