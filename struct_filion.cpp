/* Implement structural descriptors for the machine-learning technique as described in arXiv:2105.05921v1 */

#include "struct_filion.h"
#include "defs.h"
#include "pbc.h"

#include <cmath>

int NRadial=0;
int NAngular=0;
int NTot = 0;

void eval_struct_filion(){

    std::cout << "EVAL STRUCT FILION " << std::endl; 

    init_descriptors();

    int NAngular_loc = NAngular/(lmax-lmin);
    double struct_filion_save_thermal [2*(2*lmax+1)*lmax*NAngular_loc+lmax*NAngular_loc];
    //double struct_filion_save_inherent [2*(2*lmax+1)*lmax*NAngular_loc+lmax*NAngular_loc];
    int iType,jType;
    for (int s=0; s<NS;s++) { // loop over structures

        for (int i=0; i<N;i++) { // loop over particles
            //iType = type_data[i+s*N];

            for (int l=0; l< 2*(2*lmax+1)*lmax*NAngular_loc+lmax*NAngular_loc; l++) {
                struct_filion_save_thermal [l]=0.0;
                //struct_filion_save_inherent [l]=0.0;
            }

            for (int j=0; j<N;j++) { // loop over particle pairs
                if (j!=i) {
                    jType = type_data[j+s*N];
                    double dr=0.0, dx[dim];
                    double dr_inh=0.0, dx_inh[dim];
                    for (int d=0; d<dim;d++) {
                        dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        dx_inh[d] = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                        apply_pbc(dx[d]);
                        apply_pbc(dx_inh[d]);
                        dr += dx[d]*dx[d];
                        dr_inh += dx_inh[d]*dx_inh[d];
                    }

                    if (dr < 36.0) { // all other contributions are negligible
                        // eval radial descriptors
                        eval_radial(i+s*N,jType,dx,sqrt(dr),struct_filion_classifiers_thermal);
                        //eval_radial(i+s*N,jType,dx_inh,sqrt(dr_inh),struct_filion_classifiers_inherent);
                    }

                    if (dr < 10.0) {  // all other contributions are negligible
                        // eval angular descriptors
                        eval_angular(i+s*N,dx,sqrt(dr),struct_filion_save_thermal);
                        //eval_angular(i+s*N,dx_inh,sqrt(dr_inh),struct_filion_save_inherent);
                    }
                    
                } 
            }

            // collect angular descriptors
            collect_angular_descriptors(i+s*N, struct_filion_save_thermal, struct_filion_classifiers_thermal);
            //collect_angular_descriptors(i+s*N, struct_filion_save_inherent, struct_filion_classifiers_inherent);

        }
    }

    std::cout << "EVAL STRUCT FILION 2 " << std::endl; 

    normalize_cg_descriptors(0,struct_filion_classifiers_thermal);
    //normalize_cg_descriptors(1,struct_filion_classifiers_inherent);

    // print descriptors
    QString pathPred = QString::fromStdString(folderOut);
    pathPred.append("struct_filion_thermal.csv");
    QFile outfilePred3(pathPred);   // input file with xyz
    outfilePred3.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred3(&outfilePred3);
    //write header
    outPred3 << "TYPE" << ",";
    for (int cg=0; cg <= CG_NMAX; cg++) {
        for (int type=0; type<NTYPE; type++) {
            for (int k=0;k<NRadial;k++) {
                outPred3 << "CG" << cg << ":R" << struct_filion_descriptor_list[k][0] << ":" << struct_filion_descriptor_list[k][1] << ":T" << type << ",";
            }
        }
        for (int k=0;k<NAngular;k++) {
            outPred3 << "CG" << cg << ":A" << struct_filion_descriptor_list[NRadial+k][0] << ":" << struct_filion_descriptor_list[NRadial+k][1] << ":" << struct_filion_descriptor_list[NRadial+k][2];
            if (cg != CG_NMAX || k!= NAngular-1) outPred3 << ",";
        }
    }
    outPred3 << "\n";
    // write body
    for (int i=0; i<N*NS; i++) {
        outPred3 << type_data[i] << ",";
        for (int cg=0; cg <= CG_NMAX; cg++) {
            for (int k=0;k<NTot;k++) {
                outPred3 << struct_filion_classifiers_thermal[i][NTot*cg+k];
                if (cg != CG_NMAX || k!= NTot-1) outPred3 << ",";
            }
        }
        outPred3 << "\n";
    }
    outfilePred3.close();

    /*QString pathPred1 = QString::fromStdString(folderOut);
    pathPred1.append("struct_filion_inherent.csv");
    QFile outfilePred1(pathPred1);   // input file with xyz
    outfilePred1.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred1(&outfilePred1);
    //write header
    outPred1 << "TYPE" << ",";
    for (int cg=0; cg <= CG_NMAX; cg++) {
        for (int type=0; type<NTYPE; type++) {
            for (int k=0;k<NRadial;k++) {
                outPred1 << "INHCG" << cg << ":R" << struct_filion_descriptor_list[k][0] << ":" << struct_filion_descriptor_list[k][1] << ":T" << type << ",";
            }
        }
        for (int k=0;k<NAngular;k++) {
            outPred1 << "INHCG" << cg << ":A" << struct_filion_descriptor_list[NRadial+k][0] << ":" << struct_filion_descriptor_list[NRadial+k][1] << ":" << struct_filion_descriptor_list[NRadial+k][2];
            if (cg != CG_NMAX || k!= NAngular-1) outPred1 << ",";
        }
    }
    outPred1 << "\n";
    // write body
    for (int i=0; i<N*NS; i++) {
        outPred1 << type_data[i] << ",";
        for (int cg=0; cg <= CG_NMAX; cg++) {
            for (int k=0;k<NTot;k++) {
                outPred1 << struct_filion_classifiers_inherent[i][NTot*cg+k];
                if (cg != CG_NMAX || k!= NTot-1) outPred1 << ",";
            }
        }
        outPred1 << "\n";
    }
    outfilePred1.close();*/

    // free memory
    free_dmatrix(struct_filion_classifiers_thermal,0,N*NS-1,0,NTot*(CG_NMAX+1)-1);
    //free_dmatrix(struct_filion_classifiers_inherent,0,N*NS-1,0,NTot*(CG_NMAX+1)-1);

}


void init_descriptors(){

    // read descriptor file
    std::cout << "EVAL STRUCT FILION: READ DESCRIPTORS " << std::endl; 
    const QRegExp rxString(QLatin1String("[^A-Za-z0-9./_-]+"));
    QString path = QString::fromStdString(folderOut);
    path.append(QString::fromStdString(struct_filion_file));
    QFile infile(path);   
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);

    // read number of radial descriptors
    int line_counter=0;
    int NRadial_loc=0,NAngular_loc=0;
    while (!infileStream.atEnd()){
        QString line = infileStream.readLine();
        if( line.at( 0 ) == '#' ) continue;
        switch (line_counter) {
            case 0 : {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                NRadial = parts[0].toDouble();
                break;
            }
            case 1 : {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                struct_filion_descriptor_list[NRadial_loc][0] = parts[0].toDouble();
                struct_filion_descriptor_list[NRadial_loc][1] = 1.0/(2.0*parts[1].toDouble()*parts[1].toDouble());
                NRadial_loc++;
                break;
            }
            case 2 : {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                NAngular = parts[0].toDouble();
                lmin = parts[1].toInt();
                lmax = parts[2].toInt();
                break;
            }
            case 3 : {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                for (int l=lmin; l<lmax;l++) {
                    struct_filion_descriptor_list[NRadial+NAngular_loc*(lmax-lmin)+l-lmin][0] = parts[0].toDouble();
                    struct_filion_descriptor_list[NRadial+NAngular_loc*(lmax-lmin)+l-lmin][1] = 1.0/(2.0*parts[1].toDouble()*parts[1].toDouble());
                    struct_filion_descriptor_list[NRadial+NAngular_loc*(lmax-lmin)+l-lmin][2] = l;
                }
                NAngular_loc++;
                break;
            }
        }
        if (line_counter == 1) {
            // jump back if there are still more radial variables to read
            if (NRadial_loc < NRadial) line_counter--;
        }
        if (line_counter == 3) {
            // jump back if there are still more radial variables to read
            if (NAngular_loc < NAngular) line_counter--;
        }
        line_counter++;
    }

    NAngular = NAngular * (lmax-lmin);
    NTot = NTYPE*NRadial + NAngular;


    // initialize
    struct_filion_classifiers_thermal = dmatrix(0,N*NS-1,0,NTot*(CG_NMAX+1)-1);
    //struct_filion_classifiers_inherent = dmatrix(0,N*NS-1,0,NTot*(CG_NMAX+1)-1);
    for (int i=0; i<N*NS;i++) {
        for (int j = 0; j < NTot*(CG_NMAX+1); j++) {
            struct_filion_classifiers_thermal[i][j] = 0.0;
            //struct_filion_classifiers_inherent[i][j] = 0.0;
        }
    }
}

void eval_radial(int i,int jType,double * dx, double dr, double ** out){
    for (int k=0;k<NRadial;k++) {
        double dist = dr-struct_filion_descriptor_list[k][0];
        out[i][jType*NRadial+k] +=  exp( - dist * dist * struct_filion_descriptor_list[k][1] );
    }
}

void eval_angular(int i,double * dx, double dr, double * out){
    double psi[2*(2*lmax+1)*lmax];
    if (dim==2) { // just calculate the angles in 2d
        double thetaij= (dx[1] > 0.0) ? acos(dx[0]/dr) : 2*M_PI-acos(dx[0]/dr);

        for (int l=lmin; l<lmax;l++) {
            psi[2*l] = cos(l*thetaij);
            psi[2*l+1] = sin(l*thetaij);
        }

        for (int k=0;k<NAngular;k++) {
            int ind = NRadial + k;
            double dist = dr-struct_filion_descriptor_list[ind][0];
            double factor = exp( - dist * dist * struct_filion_descriptor_list[ind][1] );
            int l =  struct_filion_descriptor_list[ind][2] + 0.01;
            //std::cout << l << " " << factor*psi[2*l] << " " << psi[2*l] << " " << factor << std::endl;
            out[2*k] +=  factor*psi[2*l];
            out[2*k+1] +=  factor*psi[2*l+1];
            out[2*NAngular+k] +=  factor;
        }
    } else {
        double thetaij = acos(dx[2]/dr);
        double phiij= atan(dx[1]/dx[0]);
        if (dx[0] < 0) 
            if (dx[1]<0) phiij=phiij - M_PI;
            else phiij=phiij + M_PI;
        if (dx[0])

        // prepare spherical harmonics
        for (int l=lmin; l<lmax;l++) {
            for (int m=0; m<=l;m++) {
                double leg = sph_legendre (l, m, thetaij );
                psi[2*(2*lmax+1)*l+2*(m+lmax)] = leg * cos(m*phiij);
                psi[2*(2*lmax+1)*l+2*(m+lmax)+1] = leg * sin(m*phiij);
                //std::cout << l << " "<< m << " "<< psi[2*(2*lmax+1)*l+2*(m+lmax)] << " "<< psi[2*(2*lmax+1)*l+2*(m+lmax)+1] << std::endl;
            }
            for (int m=-l; m<0;m++) {
                double sign = 1.0;
                if ( (-m)%2==1) sign=-1.0;
                psi[2*(2*lmax+1)*l+2*(m+lmax)] = sign*psi[2*(2*lmax+1)*l+2*(-m+lmax)+1];
                psi[2*(2*lmax+1)*l+2*(m+lmax)+1] = sign*psi[2*(2*lmax+1)*l+2*(-m+lmax)];
                //std::cout << l << " "<< m << " "<< psi[2*(2*lmax+1)*l+2*(m+lmax)] << " "<< psi[2*(2*lmax+1)*l+2*(m+lmax)+1] << std::endl;
            }
        }

        for (int k=0;k<NAngular;k++) {
            int ind = NRadial + k;
            double dist = dr-struct_filion_descriptor_list[ind][0];
            double factor = exp( - dist * dist * struct_filion_descriptor_list[ind][1] );
            int l =  struct_filion_descriptor_list[ind][2] + 0.01;
            //std::cout << l << " " << factor*psi[2*l] << " " << psi[2*l] << " " << factor << std::endl;
            for (int m=-l; m<=l;m++) {
                out[2*k*(2*lmax+1)+2*(m+lmax)] +=  factor*psi[2*(2*lmax+1)*l+2*(m+lmax)];
                out[2*k*(2*lmax+1)+2*(m+lmax)+1] +=  factor*psi[2*(2*lmax+1)*l+2*(m+lmax)+1];
            }
            out[2*NAngular*(2*lmax+1)+k] +=  factor;
        }
        
    }

}

void collect_angular_descriptors(int i, double * struct_filion_save, double ** struct_filion_classifiers){

    if (dim==2) {
        for (int k=0;k<NAngular;k++) {
            double real_val = struct_filion_save[2*k];
            double imag_val = struct_filion_save[2*k+1];
            struct_filion_classifiers[i][NTYPE*NRadial+k] =  sqrt(real_val*real_val + imag_val*imag_val)/struct_filion_save[2*NAngular+k];
        }
    } else {
        for (int k=0;k<NAngular;k++) {
            int ind = NRadial + k;
            int l =  struct_filion_descriptor_list[ind][2] + 0.01;
            for (int m=-l; m<=l;m++) {
                double real_val = struct_filion_save[2*k*(2*lmax+1)+2*(m+lmax)];
                double imag_val = struct_filion_save[2*k*(2*lmax+1)+2*(m+lmax)+1];
                struct_filion_classifiers[i][NTYPE*NRadial+k] +=  (real_val*real_val + imag_val*imag_val)/struct_filion_save[2*NAngular*(2*lmax+1)+k];
            }
            struct_filion_classifiers[i][NTYPE*NRadial+k] = sqrt(4.0*M_PI/(2*l+1)*struct_filion_classifiers[i][NTYPE*NRadial+k]);
        }

    }

}

void normalize_cg_descriptors(int struct_filion_mode, double ** struct_filion_classifiers){

    // then calculate coarse-grained levels
    for (int cg=1; cg <= CG_NMAX; cg++) {

        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) { // loop over particles
                double C = 0.0;

                for (int j=0; j<N;j++) { // loop over particle pairs
                    double dr=0.0, dx[dim];
                    for (int d=0; d<dim;d++) {
                        if (struct_filion_mode == 1) dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        else dx[d] = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                        apply_pbc(dx[d]);
                        dr += dx[d]*dx[d];
                    }
                    dr = sqrt(dr);
                    double fac = exp(-dr/CG_RC);
                    //std::cout << fac << std::endl;
                    C+=fac;
                    for (int k=0; k<NTot;k++) {
                        struct_filion_classifiers[i+s*N][NTot*cg+k] += fac*struct_filion_classifiers[j+s*N][NTot*(cg-1)+k];
                    }
                }
                for (int k=0; k<NTot;k++) {
                    struct_filion_classifiers[i+s*N][NTot*cg+k] /= C;
                }
            }
        }
        //std::cout << struct_filion_classifiers[0][NTot*cg+10] << " " << struct_filion_classifiers[0][NTot*(cg-1)+10] << std::endl; 
    }

    

        std::cout << "EVAL STRUCT FILION 3 " << std::endl; 

    // normalize descriptors to have mean zero and unit variance
    for (int k=0; k<NTot*(CG_NMAX+1);k++) {
        double mean = 0.0;
        double var = 0.0;
        //std::cout << k << std::endl;
        for (int i=0; i<N*NS;i++) {
            mean += struct_filion_classifiers[i][k];
            var += struct_filion_classifiers[i][k]*struct_filion_classifiers[i][k];
        }
        mean/=N*NS;
        var=sqrt(var/(N*NS)- mean*mean);

        //std::cout << mean << " " << var << std::endl;

        for (int i=0; i<N*NS;i++) {
            struct_filion_classifiers[i][k] -= mean;
            if (var > 0.000001) struct_filion_classifiers[i][k] /= var;
        }

        // test
        /*mean = 0.0;
        var = 0.0;
        for (int i=0; i<N*NS;i++) {
            mean += struct_filion_classifiers[i][k];
            var += struct_filion_classifiers[i][k]*struct_filion_classifiers[i][k];
        }
        mean/=N*NS;
        var=sqrt(var/(N*NS)- mean*mean);
        std::cout << mean << " " << var << std::endl;*/
    }

       std::cout << "EVAL STRUCT FILION 4 " << std::endl; 

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
                for (int d=0; d<dim;d++) {
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

    QString pathOrig = QString::fromStdString(folderOut);
    QString pathPred = pathOrig;

    pathPred = pathOrig;
    pathPred.append("struct_filion_phys.csv");
    QFile outfilePred2(pathPred);   // input file with xyz
    outfilePred2.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred2(&outfilePred2);
    //write header
    outPred2 << "TYPE" << ",";
    for (int k=0; k<NStructTotal; k++) for (int c=0; c<NCG; c++) outPred2 << QString::fromStdString(StructNames[k]) << "CGP" << c << ",";
    for (int k=0; k<2*Nother; k++) for (int c=1; c<NCG; c++) outPred2 << QString::fromStdString(StructNames[struct_base_flag+2+k]) << "CGSTD" << c << ",";
    outPred2 << "\n";
    // write body
    for (int i=0; i<N*NS; i++) {
        outPred2 << type_data[i] << ",";
        for (int k=0; k<NStructTotal*NCG; k++) outPred2 << struct_local_norm[i][k] << ",";
        for (int k=0; k<2*Nother*(NCG-1); k++) outPred2 << struct_std_norm[i][k] << ",";
        outPred2 << "\n";
    }
    outfilePred2.close();
    
    free_dmatrix(struct_local_norm,0,N*NS-1,0,NStructTotal*NCG-1);
    free_dmatrix(struct_std_norm,0,N*NS-1,0,2*Nother*(NCG-1)-1);
    
}

void write_descriptors_csv_dyn(){

    QString pathOrig = QString::fromStdString(folderOut);
    QString pathPred = pathOrig;
    pathPred.append("struct_filion_labels.csv");
    QFile outfilePred(pathPred);   // input file with xyz
    outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred(&outfilePred);
    //write header
    outPred << "TYPE" << ",";
    for (int k=0; k<NDynTotal; k++) for (int t=1; t<NT; t++) outPred << QString::fromStdString(DynNames[k]) << t << ",";
    outPred << "\n";
    // write body
    for (int i=0; i<N*NS; i++) {
        outPred << type_data[i] << ",";
        for (int k=0; k<NDynTotal; k++) for (int t=1; t<NT; t++) outPred << dyn_avg_save[i][k*NT+t] << ",";
        outPred << "\n";
    }
    outfilePred.close();
    
}
