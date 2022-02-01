/* Implement structural descriptors for the machine-learning technique as described in arXiv:2105.05921v1 */

#include "struct_filion.h"
#include "defs.h"
#include "pbc.h"

int NRadial=0;
int NAngular=0;
int NTot = 0;
double psi[2*LMAX];

void eval_struct_filion(){

    std::cout << "EVAL STRUCT FILION " << std::endl; 

    // initialize
    for (int i=0; i<N*NS;i++) {
        for (int j = 0; j < 1000; j++) {
            struct_filion_classifiers[i][j] = 0.0;
        }
    }

    init_descriptors();

    int iType,jType;
    for (int s=0; s<NS;s++) { // loop over structures

        for (int i=0; i<N;i++) { // loop over particles
            //iType = type_data[i+s*N];

            for (int j=0; j<N;j++) { // loop over particle pairs
                if (j!=i) {
                    jType = type_data[j+s*N];
                    double dr=0.0, dx[dim];
                    for (int d=0; d<dim;d++) {
                        dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        apply_pbc(dx[d]);
                        dr += dx[d]*dx[d];
                    }

                    if (dr < 36.0) { // all other contributions are negligible
                        // eval radial descriptors
                        eval_radial(i+s*N,jType,dx,sqrt(dr));
                    }

                    if (dr < 10.0) {  // all other contributions are negligible
                        // eval angular descriptors
                        eval_angular(i+s*N,dx,sqrt(dr));
                    }
                    
                } 
            }

        }
    }

    std::cout << "EVAL STRUCT FILION 2 " << std::endl; 

    normalize_cg_descriptors();

    if (struct_filion_mode==1) {
        read_eval_struct_filion_ml();
    }

}


void init_descriptors(){

    // 28 radial descriptors in (0.6, 2.0]
    double delta = 0.05;
    for (int i=0; i<28;i++) {
        struct_filion_descriptor_list[NRadial][0] = 0.5 + (i+1)*delta;
        struct_filion_descriptor_list[NRadial][1] = 1.0/(2*delta*delta);
        NRadial++;
    }
    // 30 radial descriptors in (2.0, 5.0]
    delta = 0.1;
    for (int i=0; i<30;i++) {
        struct_filion_descriptor_list[NRadial][0] = 2.0 + (i+1)*delta;
        struct_filion_descriptor_list[NRadial][1] = 1.0/(2*delta*delta);
        NRadial++;
    }

    // angular descriptors
    delta = 0.1;
    for (int i=0; i<16;i++) {
        for (int l=1; l<LMAX;l++) {
            struct_filion_descriptor_list[NRadial+NAngular][0] = 1.0 + (i+1)*delta;
            struct_filion_descriptor_list[NRadial+NAngular][1] = 1.0/(2*delta*delta);
            struct_filion_descriptor_list[NRadial+NAngular][2] = l;
            NAngular++;
        }
    }

    NTot = NTYPE*NRadial + NAngular;
}

void eval_radial(int i,int jType,double * dx, double dr){
    for (int k=0;k<NRadial;k++) {
        double dist = dr-struct_filion_descriptor_list[k][0];
        struct_filion_classifiers[i][jType*NRadial+k] +=  exp( - dist * dist * struct_filion_descriptor_list[k][1] );
    }
}

void eval_angular(int i,double * dx, double dr){

    double thetaij= (dx[1] > 0.0) ? acos(dx[0]/dr) : 2*M_PI-acos(dx[0]/dr);

    for (int l=1; l<LMAX;l++) {
        psi[2*l] = cos(l*thetaij);
        psi[2*l+1] = sin(l*thetaij);
    }

    for (int k=0;k<NAngular;k++) {
        int ind = NRadial + k;
        double dist = dr-struct_filion_descriptor_list[ind][0];
        double factor = exp( - dist * dist * struct_filion_descriptor_list[ind][1] );
        int l =  struct_filion_descriptor_list[ind][2] + 0.01;
        //std::cout << l << " " << factor*psi[2*l] << " " << psi[2*l] << " " << factor << std::endl;
        struct_filion_classifiers[i][NTYPE*NRadial+k] +=  factor*psi[2*l];
        struct_filion_classifiers[i][NTYPE*NRadial+NAngular+k] +=  factor*psi[2*l+1];
        struct_filion_classifiers[i][NTYPE*NRadial+2*NAngular+k] +=  factor;
    }

}

void normalize_cg_descriptors(){

    // first calculate and normalize angular descriptors
    for (int i=0; i<N*NS;i++) {
        for (int k=0;k<NAngular;k++) {
            double real_val = struct_filion_classifiers[i][NTYPE*NRadial+k];
            double imag_val = struct_filion_classifiers[i][NTYPE*NRadial+NAngular+k];
            struct_filion_classifiers[i][NTYPE*NRadial+k] =  sqrt(real_val*real_val + imag_val*imag_val)/struct_filion_classifiers[i][NTYPE*NRadial+2*NAngular+k];
        }
        for (int k=0;k<2*NAngular;k++) {
            struct_filion_classifiers[i][NTot+k] = 0.0;
        }
    }

    // then calculate coarse-grained levels
    for (int cg=1; cg <= CG_NMAX; cg++) {

        for (int s=0; s<NS;s++) { // loop over structures
            for (int i=0; i<N;i++) { // loop over particles
                double C = 0.0;

                for (int j=0; j<N;j++) { // loop over particle pairs
                    double dr=0.0, dx[dim];
                    for (int d=0; d<dim;d++) {
                        dx[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        apply_pbc(dx[d]);
                        dr += dx[d]*dx[d];
                    }
                    dr = sqrt(dr);
                    double fac = exp(-dr/CG_RC);
                    C+=fac;
                    for (int k=0; k<NTot;k++) {
                        struct_filion_classifiers[i+s*N][NTot*cg+k] += fac*struct_filion_classifiers[i+s*N][NTot*(cg-1)+k];
                    }
                }

                for (int k=0; k<NTot;k++) {
                    struct_filion_classifiers[i+s*N][NTot*cg+k] /= C;
                }
            }
        }
    }

        std::cout << "EVAL STRUCT FILION 3 " << std::endl; 

    // normalize descriptors to have mean zero and unit variance
    for (int k=0; k<NTot*(CG_NMAX+1);k++) {
        double mean = 0.0;
        double var = 0.0;
        for (int i=0; i<N*NS;i++) {
            mean += struct_filion_classifiers[i][k];
            var += struct_filion_classifiers[i][k]*struct_filion_classifiers[i][k];
        }
        mean/=N*NS;
        var=sqrt(var/(N*NS)- mean*mean);

        for (int i=0; i<N*NS;i++) {
            struct_filion_classifiers[i][k] -= mean;
            struct_filion_classifiers[i][k] /= var;
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
                

void write_descriptors_csv(){

    std::cout << "test1" << std::endl;

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
            struct_local_norm[i][k] /= var;
        }
    }

        std::cout << "test2" << std::endl;

    int t=35;
    QString pathOrig = QString::fromStdString(folderOut);
    QString pathPred = pathOrig;
    pathPred.append("struct_filion_batch.csv");
    QFile outfilePred(pathPred);   // input file with xyz
    outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred(&outfilePred);
    //write header
    for (int k=0; k<NDyn; k++) outPred << QString::fromStdString(DynNames[k]) << ",";
    if (struct_base_flag>=0) { // include also other structural descriptors
        for (int k=0; k<NStructTotal; k++) for (int c=0; c<NCG; c++) outPred << QString::fromStdString(StructNames[k]) << c << ",";
    }
    for (int cg=0; cg <= CG_NMAX; cg++) {
        for (int type=0; type<NTYPE; type++) {
            for (int k=0;k<NRadial;k++) {
                outPred << "CG" << cg << ":R" << struct_filion_descriptor_list[k][0] << ":" << struct_filion_descriptor_list[k][1] << ":T" << type << ",";
            }
        }
        for (int k=0;k<NAngular;k++) {
            outPred << "CG" << cg << ":A" << struct_filion_descriptor_list[NRadial+k][0] << ":" << struct_filion_descriptor_list[NRadial+k][1] << ":" << struct_filion_descriptor_list[NRadial+k][2];
            if (cg != CG_NMAX || k!= NAngular-1) outPred << ",";
        }
    }
    outPred << "\n";
    // write body
    for (int i=0; i<N*NS; i++) {
        for (int k=0; k<NDyn; k++) outPred << dyn_avg_save[i][k*NT+t] << ",";
        for (int k=0; k<NStructTotal*NCG; k++) outPred << struct_local_norm[i][k] << ",";
        for (int cg=0; cg <= CG_NMAX; cg++) {
            for (int k=0;k<NTot;k++) {
                outPred << struct_filion_classifiers[i][NTot*cg+k];
                if (cg != CG_NMAX || k!= NTot-1) outPred << ",";
            }
        }
        outPred << "\n";
    }
    outfilePred.close();
    free_dmatrix(struct_local_norm,0,N*NS-1,0,NStructTotal*NCG-1);
    
            std::cout << "test3" << std::endl;
}


void read_eval_struct_filion_ml(){

}