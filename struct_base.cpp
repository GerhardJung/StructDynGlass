
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
                double dr_inherent = 0.0, dx_inherent;
                for (int d=0; d<NDim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                    dx_inherent = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                    apply_pbc(dx_inherent);
                    dr_inherent += dx_inherent*dx_inherent;
                }

                // add to gr histogram
                int bin = sqrt(dr)/boxL*2.0*NHistoGr;
                double sigma = determine_sigma(i, j);
                int bin_tot = sqrt(dr)/boxL*2.0*NHistoGr/sigma;
                if (bin > 0 && bin < NHistoGr) {
                    struct_base_gr[iType*NTYPE+jType][bin]+=1.0; // gr per type
                    struct_base_gr[NTYPE*NTYPE][bin_tot]+=1.0; // total gr
                }

                // determine neighbors
                if (dr_inherent < rcut2 && i != j ) {
                    neighbors[i+s*N][ncount] = j;
                    ncount ++;
                }

                // calculate epot
                if (i!=j) struct_local[NCG*(struct_base_flag+2)][i+s*N] += 0.5*calc_epot(i+s*N, j+s*N, dr);
                if (i!=j) struct_local[NCG*(struct_base_flag+3)][i+s*N] += 0.5*calc_epot(i+s*N, j+s*N, dr_inherent);
                
            } 

        }
    }


    // rescale and print gr
    rescale_print_gr();

    // calc psi_6/phi_5/theta_tanaka
    std::cout << "EVAL STRUCT BASE: CALC PSI " << std::endl; 
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
    double dx[NDim], dr;
    double dx_inherent[NDim], dr_inherent;
    double * save_psi = dvector(0,4*(2*(lmax)+1)*lmax-1);
    
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles
            for (int i=0; i< 4*(2*(lmax)+1)*lmax; i++) save_psi[i]=0.0;
            int n0 = 0;
            while (neighbors[i+s*N][n0] != -1 ) {
                //std::cout << n0 << " " << neighbors[i+s*N][n0] << std::endl;
                //if(n0>17) std::cout << neighbors[i+s*N][n0] << std::endl;
                dr = 0.0, dr_inherent=0.0;
                for (int d=0; d<NDim;d++) {
                    dx[d] = xyz_data[i+s*N][d] - xyz_data[neighbors[i+s*N][n0]+s*N][d];
                    apply_pbc(dx[d]);
                    dr += dx[d]*dx[d];
                    dx_inherent[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[neighbors[i+s*N][n0]+s*N][d];
                    apply_pbc(dx_inherent[d]);
                    dr_inherent += dx_inherent[d]*dx_inherent[d];
                }
                dr = sqrt(dr);
                dr_inherent = sqrt(dr_inherent);

                if (NDim==2) {
                    double thetaij= (dx[1] > 0) ? acos(dx[0]/dr) : 2*M_PI-acos(dx[0]/dr);
                    double thetaij_inherent= (dx_inherent[1] > 0) ? acos(dx_inherent[0]/dr_inherent) : 2*M_PI-acos(dx_inherent[0]/dr_inherent);
                    //std::cout << n0 << " " << theta << " " << dx[0]/sqrt(dr) << std::endl;

                    // calc psi
                    for (int psi_order=lmin; psi_order < lmax; psi_order++ ) {
                        save_psi[4*(psi_order)] += cos((psi_order)*thetaij);
                        save_psi[4*(psi_order)+1] += sin((psi_order)*thetaij);
                        save_psi[4*(psi_order)+2] += cos((psi_order)*thetaij_inherent);
                        save_psi[4*(psi_order)+3] += sin((psi_order)*thetaij_inherent);
                    }

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
                            for (int d=0; d<NDim;d++) {
                                dx[d] = xyz_data[i+s*N][d] - xyz_data[neighbors[i+s*N][n1]+s*N][d];
                                apply_pbc(dx[d]);
                                dr += dx[d]*dx[d];
                                dx_inherent[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[neighbors[i+s*N][n1]+s*N][d];
                                apply_pbc(dx_inherent[d]);
                                dr_inherent += dx_inherent[d]*dx_inherent[d];
                            }

                            // thermal
                            double thetaik= (dx[1] > 0) ? acos(dx[0]/sqrt(dr)) : 2*M_PI-acos(dx[0]/sqrt(dr));
                            double thetai = abs(thetaij - thetaik); 
                            if (thetai > M_PI) thetai = 2*M_PI -  thetai;
                            double sigma_ij = determine_sigma(i+s*N, j+s*N);    
                            double sigma_ik = determine_sigma(i+s*N, k+s*N);   
                            double sigma_kj = determine_sigma(k+s*N, j+s*N); 
                            double theta2 = acos ( (sigma_ij*sigma_ij + sigma_ik*sigma_ik - sigma_kj*sigma_kj) / (2.0*sigma_ij*sigma_ik)  );     
                            struct_local[NCG*(struct_base_flag+4)][i+s*N] +=  abs( thetai - theta2);    

                            // inherent
                            double thetaik_inherent= (dx_inherent[1] > 0) ? acos(dx_inherent[0]/sqrt(dr_inherent)) : 2*M_PI-acos(dx_inherent[0]/sqrt(dr_inherent));
                            double thetai_inherent = abs(thetaij_inherent - thetaik_inherent); 
                            if (thetai_inherent > M_PI) thetai_inherent = 2*M_PI -  thetai_inherent;  
                            struct_local[NCG*(struct_base_flag+5)][i+s*N] +=  abs( thetai_inherent - theta2);  
                            
                        }
                        n1++;
                    }

                } else {
                    calc_bond_order_3d(dx,dr,save_psi);
                    calc_bond_order_3d(dx_inherent,dr_inherent,&save_psi[2*(2*(lmax)+1)*lmax]);
                }

                n0 ++;
            }

            if (NDim==2) {
                // normalize theta tanaka
                struct_local[NCG*(struct_base_flag+4)][i+s*N] /= n0*2.0;
                struct_local[NCG*(struct_base_flag+5)][i+s*N] /= n0*2.0;
                // calculate psis of different order
                for (int psi_order=2*lmin; psi_order < 2*lmax; psi_order++ ) { // including thermal and inherent states
                    struct_local[NCG*(struct_base_flag+psi_order-2*lmin+6)][i+s*N] = sqrt(save_psi[2*psi_order]*save_psi[2*psi_order] + save_psi[2*psi_order+1]*save_psi[2*psi_order+1])/((double) n0);
                }
            } else {

                // theta tanaka not implemented for 3D
                struct_local[NCG*(struct_base_flag+4)][i+s*N] = 0.0;
                struct_local[NCG*(struct_base_flag+5)][i+s*N] = 0.0;
                for (int l=lmin; l<lmax; l++) {
                    for (int m=-l; m<=l;m++) {
                        double real_val = save_psi[2*l*(2*lmax+1)+2*(m+lmax)];
                        double imag_val = save_psi[2*l*(2*lmax+1)+2*(m+lmax)+1];
                        struct_local[NCG*(struct_base_flag+2*l-2*lmin+6)][i+s*N] +=  (real_val*real_val + imag_val*imag_val);
                    }
                    struct_local[NCG*(struct_base_flag+2*l-2*lmin+6)][i+s*N] = sqrt(4.0*M_PI/(2*l+1)*struct_local[NCG*(struct_base_flag+2*l-2*lmin+6)][i+s*N]);
                    for (int m=-l; m<=l;m++) { //inherent states
                        double real_val = save_psi[2*(2*lmax+1)*lmax+2*l*(2*lmax+1)+2*(m+lmax)];
                        double imag_val = save_psi[2*(2*lmax+1)*lmax+2*l*(2*lmax+1)+2*(m+lmax)+1];
                        struct_local[NCG*(struct_base_flag+2*l-2*lmin+7)][i+s*N] +=  (real_val*real_val + imag_val*imag_val);
                    }
                    struct_local[NCG*(struct_base_flag+2*l-2*lmin+7)][i+s*N] = sqrt(4.0*M_PI/(2*l+1)*struct_local[NCG*(struct_base_flag+2*l-2*lmin+7)][i+s*N]);
                }
            }

        }
    }
    free_dvector(save_psi,0,4*(2*lmax+1)*lmax-1);
}

// iterate over all particles to calculate coarse-grained quantities and density
void eval_den_cg(){
    double mean_den[NCG];
    double mean_den_inherent[NCG];
    int Nother = (lmax-lmin+2);
    double mean_rest[2*Nother*NCG];

    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles

            for (int c=0; c<NCG; c++) {
                mean_den[c] = 0.0;
                mean_den_inherent[c] = 0.0;
                for (int k=0; k<2*Nother; k++) {
                     mean_rest[c+NCG*k] = 0.0;   
                }
            }

            for (int j=0; j<N;j++) { // loop over particle pairs
                double dr = 0.0, dx;
                double dr_inherent = 0.0, dx_inherent;
                for (int d=0; d<NDim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                    dx_inherent = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                    apply_pbc(dx_inherent);
                    dr_inherent += dx_inherent*dx_inherent;
                }
                dr = sqrt(dr);
                dr_inherent = sqrt(dr_inherent);

                for (int c=0; c<NCG; c++) {
                    double L = c;
                    L/=2.0;
                    if (L < 0.1) L = 0.1;
                    double w = exp(-dr/L);
                    mean_den[c] += w;
                    double w_inherent = exp(-dr_inherent/L);
                    mean_den_inherent[c] += w_inherent;

                    for (int k=0; k<Nother; k++) {
                        mean_rest[c+NCG*2*k] += w*struct_local[NCG*(struct_base_flag+2*k+2)][j+s*N];
                        mean_rest[c+NCG*(2*k+1)] += w_inherent*struct_local[NCG*(struct_base_flag+2*k+3)][j+s*N];
                    }
                }
            }

            //std::cout << mean_epot[0]/mean_den[0] << " " << mean_epot[1]/mean_den[1] << " " << mean_epot[2]/mean_den[2] << " " << mean_epot[3]/mean_den[3] << std::endl;

            for (int c=0; c<NCG; c++) {
                double L = c;
                L/=2.0;
                if (L < 0.1) L = 0.1;
                struct_local[NCG*(struct_base_flag)+c][i+s*N] = mean_den[c]/((L+1.0)*(L+1.0)*(L+1.0) );
                struct_local[NCG*(struct_base_flag+1)+c][i+s*N] = mean_den_inherent[c]/((L+1.0)*(L+1.0)*(L+1.0) );
                if(c>0) {
                    for (int k=0; k<lmax-lmin+2; k++) {
                        struct_local[NCG*(struct_base_flag+2*k+2)+c][i+s*N] = mean_rest[c+NCG*2*k]/mean_den[c];
                        struct_local[NCG*(struct_base_flag+2*k+3)+c][i+s*N] = mean_rest[c+NCG*(2*k+1)]/mean_den_inherent[c];
                    }
                }
            }

        }
    }

}

void write_descriptors_csv_phys(){

    // calculate additional structural descriptors
    double ** struct_local_var; 
    int Nother = (lmax-lmin+2);
    struct_local_var = dmatrix(0,N*NS-1,0,2*Nother*(NCG-1)-1);
    double mean_den[NCG];
    double mean_den_inherent[NCG];
    double mean_rest[2*Nother*NCG];

    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles

            for (int c=0; c<NCG; c++) {
                mean_den[c] = 0.0;
                mean_den_inherent[c] = 0.0;
                for (int k=0; k<2*Nother; k++) {
                     mean_rest[c+NCG*k] = 0.0;   
                }
            }
            //double mean=0.0;
            for (int j=0; j<N;j++) { // loop over particle pairs
                double dr = 0.0, dx;
                double dr_inherent = 0.0, dx_inherent;
                for (int d=0; d<NDim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                    dx_inherent = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                    apply_pbc(dx_inherent);
                    dr_inherent += dx_inherent*dx_inherent;
                }
                dr = sqrt(dr);
                dr_inherent = sqrt(dr_inherent);

                for (int c=0; c<NCG; c++) {
                    double L = c;
                    L/=2.0;
                    if (L < 0.1) L = 0.1;
                    double w = exp(-dr/L);
                    mean_den[c] += w;
                    double w_inherent = exp(-dr_inherent/L);
                    mean_den_inherent[c] += w_inherent;

                    for (int k=0; k<Nother; k++) {
                        /*if (s==0 && i==0 && c==2 && k==0 && w>0.01) {
                            std::cout << w << " " << struct_local[NCG*(struct_base_flag+2*k+2)][j+s*N] << " " << struct_local[NCG*(struct_base_flag+2*k+2)+c][i+s*N] << std::endl;
                            mean += w*struct_local[NCG*(struct_base_flag+2*k+2)][j+s*N];
                        }*/
                        mean_rest[c+NCG*2*k] += w*(struct_local[NCG*(struct_base_flag+2*k+2)][j+s*N]-struct_local[NCG*(struct_base_flag+2*k+2)+c][i+s*N])*(struct_local[NCG*(struct_base_flag+2*k+2)][j+s*N]-struct_local[NCG*(struct_base_flag+2*k+2)+c][i+s*N]);
                        mean_rest[c+NCG*(2*k+1)] += w_inherent*(struct_local[NCG*(struct_base_flag+2*k+3)][j+s*N]-struct_local[NCG*(struct_base_flag+2*k+3)+c][i+s*N])*(struct_local[NCG*(struct_base_flag+2*k+3)][j+s*N]-struct_local[NCG*(struct_base_flag+2*k+3)+c][i+s*N]);
                    }
                }
            }
            //if (s==0 && i==0) std::cout << mean/mean_den[2] << std::endl;

            //std::cout << mean_epot[0]/mean_den[0] << " " << mean_epot[1]/mean_den[1] << " " << mean_epot[2]/mean_den[2] << " " << mean_epot[3]/mean_den[3] << std::endl;

            for (int c=1; c<NCG; c++) {
                double L = c;
                L/=2.0;
                for (int k=0; k<Nother; k++) {
                    struct_local_var[i+s*N][(NCG-1)*(2*k)+c-1] = mean_rest[c+NCG*2*k]/mean_den[c];
                    // calculate standard deviation
                    struct_local_var[i+s*N][(NCG-1)*(2*k)+c-1] = sqrt(struct_local_var[i+s*N][(NCG-1)*(2*k)+c-1]);

                    struct_local_var[i+s*N][(NCG-1)*(2*k+1)+c-1] = mean_rest[c+NCG*(2*k+1)]/mean_den_inherent[c];
                    struct_local_var[i+s*N][(NCG-1)*(2*k+1)+c-1] = sqrt(struct_local_var[i+s*N][(NCG-1)*(2*k+1)+c-1]);
                }
            }

            

        }
    }

    // normalize physical structural descriptors to have mean zero and unit variance
    double ** struct_local_norm; 
    struct_local_norm = dmatrix(0,N*NS-1,0,NStructTotal*NCG-1);
    for (int k=0; k<NStructTotal*NCG;k++) {
        double mean = 0.0;
        double var = 0.0;
        for (int i=0; i<N*NS;i++) {
            mean += struct_local[k][i];
            var += struct_local[k][i]*struct_local[k][i];
        }
        mean/=N*NS;
        var=sqrt(var/(N*NS)- mean*mean);

        for (int i=0; i<N*NS;i++) {
            struct_local_norm[i][k] = struct_local[k][i] - mean;
            if (var > 0.00000001) struct_local_norm[i][k] /= var;
        }
    }
    double ** struct_std_norm; 
    struct_std_norm = dmatrix(0,N*NS-1,0,2*Nother*(NCG-1)-1);
    for (int k=0; k<2*Nother*(NCG-1);k++) {
        double mean = 0.0;
        double var = 0.0;
        for (int i=0; i<N*NS;i++) {
            mean += struct_local_var[i][k];
            var += struct_local_var[i][k]*struct_local_var[i][k];
        }
        mean/=N*NS;
        var=sqrt(var/(N*NS)- mean*mean);

        for (int i=0; i<N*NS;i++) {
            struct_std_norm[i][k] = struct_local_var[i][k] - mean;
            if (var > 0.00000001) struct_std_norm[i][k] /= var;
        }
    }

    for (int type=0; type<NTYPE; type++) {
        QString pathOrig = QString::fromStdString(folderOut);
        QString pathPred = pathOrig;
        pathPred.append("struct_filion_phys_type");
        pathPred.append(QString("%1.csv").arg(type));
        QFile outfilePred2(pathPred);   // input file with xyz
        outfilePred2.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred2(&outfilePred2);
        //write header
        for (int k=0; k<NStructTotal; k++) for (int c=0; c<NCG; c++) outPred2 << QString::fromStdString(StructNames[k]) << "CGP" << c << ",";
        for (int k=0; k<2*Nother; k++) for (int c=1; c<NCG; c++) outPred2 << QString::fromStdString(StructNames[struct_base_flag+2+k]) << "CGSTD" << c << ",";
        outPred2 << "\n";
        // write body
        for (int i=0; i<N*NS; i++) {
            if (type_data[i] == type) {
                for (int k=0; k<NStructTotal*NCG; k++) outPred2 << struct_local_norm[i][k] << ",";
                for (int k=0; k<2*Nother*(NCG-1); k++) outPred2 << struct_std_norm[i][k] << ",";
                outPred2 << "\n";
            }
        }
        outfilePred2.close();
    }
    
    free_dmatrix(struct_local_norm,0,N*NS-1,0,NStructTotal*NCG-1);
    free_dmatrix(struct_std_norm,0,N*NS-1,0,2*Nother*(NCG-1)-1);
    
}


// help function
double determine_sigma(int i, int j) {
    if(model=="KA2-2D") {
        int iType = type_data[i];
        int jType = type_data[j];
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.88;
        if (iType==2 && jType == 2) return 0.94;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 0.8;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.9;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 0.8;
    }
    if (model=="KA2" || model=="KA") {
        int iType = type_data[i];
        int jType = type_data[j];
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.88;
        if (iType==2 && jType == 2) return 0.94;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 0.8;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.9;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 0.84;
    }
    if (model=="POLY") {
        double iRadius = dia_data[i];
        double jRadius = dia_data[j];
        return (iRadius+jRadius)/2.0*(1.0-0.2*abs(iRadius-jRadius));
    }
    return 0.0;
}
double determine_epsilon(int iType, int jType) {
    if(model=="KA2-2D") {
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.5;
        if (iType==2 && jType == 2) return 0.75;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 1.5;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.75;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 1.5;
    } 
    if(model=="KA2" || model=="KA") {
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.5;
        if (iType==2 && jType == 2) return 0.75;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 1.5;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 1.25;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 1.0;
    } 
    if(model=="POLY") {
        return 1.0;
    } 
    return 0.0;
}

double calc_epot(int i, int j, double dist2) {
    double sigma = determine_sigma(i, j);
    double epsilon = determine_epsilon(type_data[i], type_data[j]);
    double sigma2=sigma*sigma;

    if(model=="KA2" || model=="KA2-2D") {
        if (dist2 < RC2LJ*sigma2 ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(C0LJ+C2LJ*rij2+C4LJ*rij4 - rij6 + rij12 ) ;
        } 
    } if(model=="KA") {
        if (dist2 < RC2LJ ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(- rij6 + rij12 ) ;
        } 
    } if(model=="POLY") {
        if (dist2 < RC2POLY*sigma2 ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(C0POLY+C2POLY*rij2+C4POLY*rij4 + rij12 ) ;
        } 
    } else return 0.0;
}

void calc_bond_order_3d(double * dx,double dr, double * save_psi){
    double thetaij = acos(dx[2]/dr);
    double phiij= atan(dx[1]/dx[0]);
    if (dx[0] < 0) 
        if (dx[1]<0) phiij=phiij - M_PI;
        else phiij=phiij + M_PI;

    // prepare spherical harmonics
    for (int l = lmin; l < lmax; l++) {
        for (int m=0; m<=l;m++) {
            double leg = sph_legendre (l, m, thetaij );
            save_psi[2*(2*lmax+1)*l+2*(m+lmax)] += leg * cos(m*phiij);
            save_psi[2*(2*lmax+1)*l+2*(m+lmax)+1] += leg * sin(m*phiij);
            //std::cout << l << " "<< m << " "<< lmax << " "<< 2*(2*lmax+1)*l+2*(m+lmax) << " "<< 4*(2*lmax+1)*lmax-1 << std::endl;
            if (save_psi[2*(2*lmax+1)*l+2*(m+lmax)]!=save_psi[2*(2*lmax+1)*l+2*(m+lmax)]) {
                //std::cout << l << " " << m << " " << dx[2]/dr << " " << phiij << " " << save_psi[2*(2*lmax+1)*l+2*(m+lmax)] << std::endl;
                exit(0);
            }
        }
        for (int m=-l; m<0;m++) {
            double sign = 1.0;
            if ( (-m)%2==1) sign=-1.0;
            save_psi[2*(2*lmax+1)*l+2*(m+lmax)] += sign*save_psi[2*(2*lmax+1)*l+2*(-m+lmax)+1];
            save_psi[2*(2*lmax+1)*l+2*(m+lmax)+1] += sign*save_psi[2*(2*lmax+1)*l+2*(-m+lmax)];
            //std::cout << l << " "<< m << " "<< psi[2*(2*lmax+1)*l+2*(m+lmax)] << " "<< psi[2*(2*lmax+1)*l+2*(m+lmax)+1] << std::endl;
        }
    }

}