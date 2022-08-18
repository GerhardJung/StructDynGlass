/* Implement structural descriptors for the machine-learning technique as described in arXiv:2105.05921v1 */

#include "struct_filion.h"
#include "struct_base.h"
#include "defs.h"
#include "pbc.h"
#include "voro++_2d.hh"
using namespace voro;

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
                    double dr=0.0, dx[NDim];
                    double dr_inh=0.0, dx_inh[NDim];
                    for (int d=0; d<NDim;d++) {
                        dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        dx_inh[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        apply_pbc(dx[d]);
                        apply_pbc(dx_inh[d]);
                        dr += dx[d]*dx[d];
                        dr_inh += dx_inh[d]*dx_inh[d];
                    }
                    dr = sqrt(dr);

                    if (dr < 6.0) { // all other contributions are negligible
                        // eval radial descriptors
                        eval_radial(i+s*N,jType,dx,dr,struct_filion_classifiers_thermal);
                        //eval_radial(i+s*N,jType,dx_inh,sqrt(dr_inh),struct_filion_classifiers_inherent);
                    }

                    if (dr < 1.5) {  // all other contributions are negligible
                        // eval angular descriptors
                        eval_angular(i+s*N,dx,dr,struct_filion_save_thermal);
                        //eval_angular(i+s*N,dx_inh,sqrt(dr_inh),struct_filion_save_inherent);
                    }
                    
                } 
            }
            //if (i<5) std::cout << struct_local_filion[NCG][i+s*N] << std::endl;

            // collect angular descriptors
            collect_angular_descriptors(i+s*N, struct_filion_save_thermal, struct_filion_classifiers_thermal);
            //collect_angular_descriptors(i+s*N, struct_filion_save_inherent, struct_filion_classifiers_inherent);

        }

    }

    std::cout << "EVAL STRUCT FILION 2 " << std::endl; 

    normalize_cg_descriptors(0,struct_filion_classifiers_thermal);
    //normalize_cg_descriptors(1,struct_filion_classifiers_inherent);

    // print descriptors
    for (int type=0; type<NTYPE; type++) {
        QString pathPred = QString::fromStdString(folderOut);
        pathPred.append("struct_filion_thermal_type");
        pathPred.append(QString("%1.csv").arg(type));
        QFile outfilePred3(pathPred);   // input file with xyz
        outfilePred3.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred3(&outfilePred3);
        //write header
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
        // length scale
        for (int cg=0; cg <= CG_NMAX; cg++) {
            for (int type=0; type<NTYPE; type++) {
                for (int k=0;k<NRadial;k++) {
                    outPred3 <<  struct_filion_descriptor_list[k][0] +cg*CG_RC << ",";
                }
            }
            for (int k=0;k<NAngular;k++) {
                outPred3 << struct_filion_descriptor_list[NRadial+k][0] +cg*CG_RC ;
                if (cg != CG_NMAX || k!= NAngular-1) outPred3 << ",";
            }
        }
        outPred3 << "\n";
        // mean and variance
        for (int k=0; k<(CG_NMAX+1)*NTot; k++) {
            outPred3 << struct_mean_var[k][2*type];
            if (k<(CG_NMAX+1)*NTot-1) outPred3 << ",";
        }
        outPred3 << "\n";
        for (int k=0; k<(CG_NMAX+1)*NTot; k++) {
            outPred3 << struct_mean_var[k][2*type+1];
            if (k<(CG_NMAX+1)*NTot-1) outPred3 << ",";
        }
        outPred3 << "\n";

        // particle data
        for (int i=0; i<N*NS; i++) {
            if (type_data[i] == type) {
                for (int cg=0; cg <= CG_NMAX; cg++) {
                    for (int k=0;k<NTot;k++) {
                        outPred3 << struct_filion_classifiers_thermal[i][NTot*cg+k];
                        if (cg != CG_NMAX || k!= NTot-1) outPred3 << ",";
                    }
                }
                outPred3 << "\n";
            }
        }
        outfilePred3.close();
    }

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

    std::cout << "EVAL STRUCT FILION 5 " << std::endl; 

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

    struct_mean_var = dmatrix(0,NTot*(CG_NMAX+1),0,2*NTYPE-1);
}

void eval_radial(int i,int jType,double * dx, double dr, double ** out){
    for (int k=0;k<NRadial;k++) {
        //double dist = dr-struct_filion_descriptor_list[k][0];
        //out[i][jType*NRadial+k] +=  exp( - dist * dist * struct_filion_descriptor_list[k][1] );
        if (dr < struct_filion_descriptor_list[k][0]) out[i][jType*NRadial+k] += 1.0;
    }
}

void eval_angular(int i,double * dx, double dr, double * out){
    double psi[2*(2*lmax+1)*lmax];
    if (NDim==2) { // just calculate the angles in 2d
        double thetaij= (dx[1] > 0.0) ? acos(dx[0]/dr) : 2*M_PI-acos(dx[0]/dr);

        for (int l=lmin; l<lmax;l++) {
            psi[2*l] = cos(l*thetaij);
            psi[2*l+1] = sin(l*thetaij);
        }

        for (int k=0;k<NAngular;k++) {
            int ind = NRadial + k;
            double dist = dr-struct_filion_descriptor_list[ind][0];
            //double factor = exp( - dist * dist * struct_filion_descriptor_list[ind][1] );
            double factor = 1.0;
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
            //double factor = exp( - dist * dist * struct_filion_descriptor_list[ind][1] );
            double factor = 1.0;
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

    if (NDim==2) {
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
            struct_filion_classifiers[i][NTYPE*NRadial+k] = sqrt(4.0*M_PI/(2.0*l+1)*struct_filion_classifiers[i][NTYPE*NRadial+k]);
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
                    double dr=0.0, dx[NDim];
                    for (int d=0; d<NDim;d++) {
                        if (struct_filion_mode == 1) dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        else dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        apply_pbc(dx[d]);
                        dr += dx[d]*dx[d];
                    }
                    //if (dr < CG_RC*CG_RC) {
                        dr = sqrt(dr);
                        double fac = exp(-dr/CG_RC);
                        //std::cout << fac << std::endl;
                        C+=fac;
                        for (int k=0; k<NTot;k++) {
                            struct_filion_classifiers[i+s*N][NTot*cg+k] += fac*struct_filion_classifiers[j+s*N][NTot*(cg-1)+k];
                        }
                    //}
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
        double mean[NTYPE] = {0.0};
        double var[NTYPE] = {0.0};
        //std::cout << k << std::endl;
        for (int i=0; i<N*NS;i++) {
            //if (k==0 && i<1000) std::cout << struct_filion_classifiers[i][k] << std::endl;
            int type = type_data[i];
            mean[type] += struct_filion_classifiers[i][k];
            var[type] += struct_filion_classifiers[i][k]*struct_filion_classifiers[i][k];
        }

        for (int type=0; type<NTYPE; type++) {
            mean[type]/=NPerType[type]*NS;
            var[type]=sqrt(var[type]/(NPerType[type]*NS)- mean[type]*mean[type]);
            struct_mean_var[k][2*type] = mean[type];
            struct_mean_var[k][2*type+1] = var[type];
        }

        //if (k==0) std::cout << mean << " " << var << std::endl;

        
        for (int i=0; i<N*NS;i++) {
            //struct_filion_classifiers[i][k] -= mean;
            //if (var > 0.000001) struct_filion_classifiers[i][k] /= var;
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
