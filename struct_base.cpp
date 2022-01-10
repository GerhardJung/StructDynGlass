
#include "struct_base.h"
#include "defs.h"
#include "dyn_bb.h"
#include "pbc.h"

void eval_struct_base(){

    std::cout << "EVAL STRUCT BASE " << std::endl; 

    // save neighbors
    int ** neighbors = imatrix(0,NS*N-1,0,N_NEIGH_MAX-1);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            for (int d=0; d<N_NEIGH_MAX;d++) {
                neighbors[i+s*N][d] = -1;
            }
        }
    }
    int iType,jType;

    for (int s=0; s<NS;s++) { // loop over structures

        for (int i=0; i<N;i++) { // loop over particles
            int ncount = 0;
            iType = type_data[i+s*N];

            for (int j=0; j<N;j++) { // loop over particle pairs
                jType = type_data[j+s*N];
                double dr = 0.0, dx;
                for (int d=0; d<dim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                }

                // add to gr histogram
                int bin = sqrt(dr)/boxL*2.0*NHistoGr;
                if (bin > 0 && bin < NHistoGr) {
                    struct_base_gr[iType*NTYPE+jType][bin]+=1.0; // gr per type
                    struct_base_gr[NTYPE*NTYPE][bin]+=1.0; // total gr
                }

                // determine neighbors
                if (dr < rcut2 && i != j ) {
                    neighbors[i+s*N][ncount] = j;
                    ncount ++;
                }

                // calculate epot
                if (i!=j) struct_local[NCG*(struct_base_flag+1)][i+s*N] += calc_epot(iType, jType, dr);
                
            } 

        }
    }


    // rescale and print gr
    rescale_print_gr();

    // calc psi_6/phi_5/theta_tanaka
    calc_psi(neighbors);

    // coarse-grain and calculate density
    eval_den_cg();

    free_imatrix(neighbors,0,NS*N-1,0,N_NEIGH_MAX-1);
}


void rescale_print_gr(){
    // rescale gr
    for (int k=0; k<NHistoGr; k++) {
        double rup = (k + 1.0) / ((double) NHistoGr ) * boxL / 2.0;
        double rlow = (k) / ((double) NHistoGr ) * boxL / 2.0;
        double dr = M_PI * (rup*rup - rlow*rlow);
        for (int i = 0; i < NTYPE; i++) {
            for (int j = 0; j < NTYPE; j++) {
                struct_base_gr[i*NTYPE+j][k] /= dr ;
                struct_base_gr[i*NTYPE+j][k] *= boxL * boxL ;
                struct_base_gr[i*NTYPE+j][k] /= (double) NPerType[i] * NPerType[j] * NS ;
            }
        }
        struct_base_gr[NTYPE*NTYPE][k] /= dr ;
        struct_base_gr[NTYPE*NTYPE][k] *= boxL * boxL ;
        struct_base_gr[NTYPE*NTYPE][k] /= (double) N * N * NS;
    }

    // print gr
    QString pathOrig = QString::fromStdString(folderOut);
    pathOrig.append("/gr.rdf");
    QFile outfileGr(pathOrig);   // input file with xyz
    outfileGr.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outGr(&outfileGr);
    outGr << "r All-All ";
    for (int i = 0; i < NTYPE; i++) {
        for (int j = 0; j < NTYPE; j++) {
            outGr << i+1 << "-" << j+1 << " ";
        }
    }
    outGr << "\n";
    for (int k=0; k<NHistoGr; k++) {
        outGr << (k + 0.5) / NHistoGr * boxL / 2.0 << " " << struct_base_gr[NTYPE*NTYPE][k] << " ";
        for (int i = 0; i < NTYPE; i++) {
            for (int j = 0; j < NTYPE; j++) {
                outGr << struct_base_gr[i*NTYPE+j][k] << " ";
            }
        }
        outGr << "\n";
    }
    outfileGr.close();

}

void calc_psi(int ** neighbors){
    double dx[dim], dr;
    double save_psi[4*N*NS] = {0};
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles
            int n0 = 0;
            while (neighbors[i+s*N][n0] != -1 ) {
                //std::cout << n0 << " " << neighbors[i+s*N][n0] << std::endl;
                //if(n0>17) std::cout << neighbors[i+s*N][n0] << std::endl;
                dr = 0.0;
                for (int d=0; d<dim;d++) {
                    dx[d] = xyz_data[i+s*N][d] - xyz_data[neighbors[i+s*N][n0]+s*N][d];
                    apply_pbc(dx[d]);
                    dr += dx[d]*dx[d];
                }
                double thetaij= (dx[1] > 0) ? acos(dx[0]/sqrt(dr)) : 2*M_PI-acos(dx[0]/sqrt(dr));
                //std::cout << n0 << " " << theta << " " << dx[0]/sqrt(dr) << std::endl;

                // calc psi
                save_psi[4*(i+s*N)] += cos(5.0*thetaij);
                save_psi[4*(i+s*N)+1] += sin(5.0*thetaij);
                save_psi[4*(i+s*N)+2] += cos(6.0*thetaij);
                save_psi[4*(i+s*N)+3] += sin(6.0*thetaij);

                // calc theta tanaka
                int n1=0;
                while (neighbors[i+s*N][n1] != -1 ) {
                    int j=neighbors[i+s*N][n0];
                    int k=neighbors[i+s*N][n1];
                    // check that n1 and n0 are also neighbors
                    int n2=0;
                    int neigh_check = 0;
                    while (neighbors[j+s*N][n2] != -1 ) {
                        if (neighbors[j+s*N][n2] == k) {
                            neigh_check=1;
                            break;
                        }
                        n2++;
                    }
                    if (neigh_check==1) { // found triangle ijk
                        dr = 0.0;
                        for (int d=0; d<dim;d++) {
                            dx[d] = xyz_data[i+s*N][d] - xyz_data[neighbors[i+s*N][n1]+s*N][d];
                            apply_pbc(dx[d]);
                            dr += dx[d]*dx[d];
                        }
                        double thetaik= (dx[1] > 0) ? acos(dx[0]/sqrt(dr)) : 2*M_PI-acos(dx[0]/sqrt(dr));
                        double thetai = abs(thetaij - thetaik); 
                        if (thetai > M_PI) thetai = 2*M_PI -  thetai;

                        double sigma_ij = determine_sigma(type_data[i+s*N], type_data[j+s*N]);    
                        double sigma_ik = determine_sigma(type_data[i+s*N], type_data[k+s*N]);   
                        double sigma_kj = determine_sigma(type_data[k+s*N], type_data[j+s*N]); 
                        double theta2 = acos ( (sigma_ij*sigma_ij + sigma_ik*sigma_ik - sigma_kj*sigma_kj) / (2.0*sigma_ij*sigma_ik)  );     
                        struct_local[NCG*(struct_base_flag+4)][i+s*N] +=  abs( thetai - theta2);    
                        //std::cout << thetai << " " << theta2 << " " << abs( thetai - theta2) << std::endl;    
                    }
                    n1++;
                }
                
                n0 ++;
            }

            for (int j = 0; j < 2; j++) {
                struct_local[NCG*(struct_base_flag+j+2)][i+s*N] = sqrt(save_psi[4*(i+s*N)+2*j]*save_psi[4*(i+s*N)+2*j] + save_psi[4*(i+s*N)+2*j+1]*save_psi[4*(i+s*N)+2*j+1])/((double) n0);
            }
            struct_local[NCG*(struct_base_flag+4)][i+s*N] /= n0*2.0;

        }
    }
}

// iterate over all particles to calculate coarse-grained quantities and density
void eval_den_cg(){
    double mean_den[NCG];
    double mean_epot[NCG];
    double mean_psi[2*NCG];
    double mean_theta_tanaka[NCG];

    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles

            for (int c=0; c<NCG; c++) {
                mean_den[c] = 0.0;
                mean_epot[c]=0.0;
                mean_psi[2*c]=0.0;
                mean_psi[2*c+1]=0.0;
                mean_theta_tanaka[c]=0.0;
            }

            for (int j=0; j<N;j++) { // loop over particle pairs
                double dr = 0.0, dx;
                for (int d=0; d<dim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                }
                dr = sqrt(dr);

                for (int c=0; c<NCG; c++) {
                    double L = c;
                    if (L < 0.1) L = 0.1;
                    double w = exp(-dr/L);
                    mean_den[c] += w;
                    mean_epot[c] += w*struct_local[NCG*(struct_base_flag+1)][j+s*N]/2.0;
                    mean_psi[2*c] += w*struct_local[NCG*(struct_base_flag+2)][j+s*N];
                    mean_psi[2*c+1] += w*struct_local[NCG*(struct_base_flag+3)][j+s*N];
                    mean_theta_tanaka[c]+= w*struct_local[NCG*(struct_base_flag+4)][j+s*N];
                }
            }

            for (int c=0; c<NCG; c++) {
                double L = c;
                if (L < 0.1) L = 0.1;
                struct_local[NCG*(struct_base_flag)+c][i+s*N] = mean_den[c]/((L+1.0)*(L+1.0)*(L+1.0) );
                struct_local[NCG*(struct_base_flag+1)+c][i+s*N] = mean_epot[c]/mean_den[c];
                struct_local[NCG*(struct_base_flag+2)+c][i+s*N] = mean_psi[2*c]/mean_den[c];
                struct_local[NCG*(struct_base_flag+3)+c][i+s*N] = mean_psi[2*c+1]/mean_den[c];
                struct_local[NCG*(struct_base_flag+4)+c][i+s*N] = mean_theta_tanaka[c]/mean_den[c];
            }

        }
    }

}


// help function
double determine_sigma(int iType, int jType) {
    if (iType==0 && jType == 0) return 1.0;
    if (iType==1 && jType == 1) return 0.88;
    if (iType==2 && jType == 2) return 0.94;
    if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 0.8;
    if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.9;
    if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 0.8;
    return 0.0;
}
double determine_epsilon(int iType, int jType) {
    if (iType==0 && jType == 0) return 1.0;
    if (iType==1 && jType == 1) return 0.5;
    if (iType==2 && jType == 2) return 0.75;
    if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 1.5;
    if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.75;
    if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 1.5;
    return 0.0;
}

double calc_epot(int iType, int jType, double dist2) {
    double sigma = determine_sigma(iType, jType);
    double epsilon = determine_epsilon(iType, jType);

    if (dist2 < RC2 ) {
        double rij2 = dist2/(sigma*sigma);
        double rij4 = rij2*rij2;
        double rij6 = 1.0/(rij4*rij2);
        double rij12 = rij6*rij6;
        return 4.0*epsilon*(C0+C2*rij2+C4*rij4 - rij6 + rij12 ) ;
    } else return 0.0;
}