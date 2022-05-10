/* Implement structural descriptors for the machine-learning technique as described in arXiv:2105.05921v1 */

#include "struct_ml.h"
#include "struct_base.h"
#include "defs.h"
#include "pbc.h"
#include "voro++_2d.hh"
using namespace voro;

#include <cmath>

void eval_struct_ml(){

    std::cout << "EVAL STRUCT ML" << std::endl; 

    int iType,jType;
    for (int s=0; s<NS;s++) { // loop over structures

        for (int i=0; i<N;i++) { // loop over particles
            //iType = type_data[i+s*N];

            for (int j=0; j<N;j++) { // loop over particle pairs
                if (j!=i) {
                    jType = type_data[j+s*N];
                    double dr=0.0, dx[dim];
                    double dr_inh=0.0, dx_inh[dim];
                    for (int d=0; d<dim;d++) {
                        dx_inh[d] = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                        apply_pbc(dx_inh[d]);
                        dr_inh += dx_inh[d]*dx_inh[d];
                    }
                    
                    // calculate epot
                    if (dr_inh < RC2LJ) struct_local_ml[(NTYPE+1)*NCG][i+s*N] += 0.5*calc_epot(i+s*N, j+s*N, dr_inh);
                } 
            }
        }

        // include voronoi
        calc_voronoi(s);
    }

    std::cout << "EVAL STRUCT ML 2 " << std::endl; 

    // coarse-grain and calculate density
    eval_den_cg_ml();

    // write physical descriptors
    write_descriptors_csv();

    std::cout << "EVAL STRUCT ML 3 " << std::endl; 

}

// calculate voronoi properties
void calc_voronoi(int s){

    container_poly_2d con(-boxL/2.0,boxL/2.0,-boxL/2.0,boxL/2.0,10,10,true,true,16);

    // Add particles to the container
	for(int i=0;i<N;i++) {
		double sigma = 0.5;
		if (type_data[s*N+i] == 1) sigma=0.44;
		if (type_data[s*N+i] == 2) sigma=0.47;
		con.put(i,xyz_inherent_data[s*N+i][0],xyz_inherent_data[s*N+i][1],sigma);
	}
	
	c_loop_all_2d vl(con);
	voronoicell_neighbor_2d c;
	int ij,q, id;double *pp;
	double cx,cy;
	vector<int> vi;
	if(vl.start()) do if(con.compute_cell(c,vl)) {
			ij=vl.ij;q=vl.q;pp=con.p[ij]+con.ps*q;
			c.centroid(cx,cy);
            id = con.id[ij][q];
            struct_local_ml[(2*(NTYPE+1))*NCG][id+s*N] = c.perimeter();
	} while(vl.inc());

}

// iterate over all particles to calculate coarse-grained quantities and density
void eval_den_cg_ml(){
    double mean_den_inherent[(NTYPE+1)*NCG];
    double mean_rest[4*(NTYPE+1)*NCG];

    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles

            for (int c=0; c<NCG; c++) {
                for (int type=0; type<(NTYPE+1); type++) mean_den_inherent[type*NCG+c] = 0.0;
                for (int k=0; k<4*(NTYPE+1); k++) {
                    mean_rest[k*NCG+c] = 0.0;   
                }
            }

            for (int j=0; j<N;j++) { // loop over particle pairs
                double dr_inherent = 0.0, dx_inherent;
                for (int d=0; d<dim;d++) {
                    dx_inherent = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                    apply_pbc(dx_inherent);
                    dr_inherent += dx_inherent*dx_inherent;
                }
                dr_inherent = sqrt(dr_inherent);

                for (int c=1; c<NCG; c++) {
                    double L = c;
                    L/=2.0;
                    //if (dr_inherent < L) {
                    double w_inherent = exp(-dr_inherent/L);
                    for (int type=0; type<NTYPE; type++) {
                        if (type==type_data[j+s*N]) {
                            mean_den_inherent[type*NCG+c] += w_inherent;
                            for (int k=0; k<2; k++) mean_rest[type*NCG+k*(NTYPE+1)*NCG+c] += w_inherent*struct_local_ml[(NTYPE+1)*NCG+k*(NTYPE+1)*NCG][j+s*N];
                        }
                    }
                    mean_den_inherent[NTYPE*NCG+c] += w_inherent;
                    for (int k=0; k<2; k++) mean_rest[NTYPE*NCG+k*(NTYPE+1)*NCG+c] += w_inherent*struct_local_ml[(NTYPE+1)*NCG+k*(NTYPE+1)*NCG][j+s*N];
                    //}
                }
            }

            //std::cout << mean_epot[0]/mean_den[0] << " " << mean_epot[1]/mean_den[1] << " " << mean_epot[2]/mean_den[2] << " " << mean_epot[3]/mean_den[3] << std::endl;
            for (int type=0; type<(NTYPE+1); type++) struct_local_ml[type*NCG][i+s*N] = 0.0;
            for (int c=1; c<NCG; c++) {
                double L = c;
                L/=2.0;
                for (int type=0; type<(NTYPE+1); type++) {
                    struct_local_ml[type*NCG+c][i+s*N] = mean_den_inherent[type*NCG+c];
                    for (int k=0; k<2; k++)  struct_local_ml[(NTYPE+1)*NCG+type*NCG+k*(NTYPE+1)*NCG+c][i+s*N] = mean_rest[type*NCG+k*(NTYPE+1)*NCG+c]/mean_den_inherent[type*NCG+c];
                }
            }

        }
    }

}                

void write_descriptors_csv(){

    // calculate additional structural descriptors
    double ** struct_local_var; 
    double mean_den_inherent[(NTYPE+1)*NCG];
    double mean_rest[4*(NTYPE+1)*NCG];

    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles

            for (int c=0; c<NCG; c++) {
                for (int type=0; type<(NTYPE+1); type++) mean_den_inherent[type*NCG+c] = 0.0;
                for (int k=0; k<4*(NTYPE+1); k++) {
                     mean_rest[k*NCG+c] = 0.0;   
                }
            }
            //double mean=0.0;
            for (int j=0; j<N;j++) { // loop over particle pairs
                double dr_inherent = 0.0, dx_inherent;
                for (int d=0; d<dim;d++) {
                    dx_inherent = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                    apply_pbc(dx_inherent);
                    dr_inherent += dx_inherent*dx_inherent;
                }
                dr_inherent = sqrt(dr_inherent);

                for (int c=1; c<NCG; c++) {
                    double L = c;
                    L/=2.0;
                    double w_inherent = exp(-dr_inherent/L);
                    for (int type=0; type<NTYPE; type++) {
                        if (type==type_data[j+s*N]) {
                            mean_den_inherent[type*NCG+c] += w_inherent;

                            double diff = (struct_local_ml[(NTYPE+1)*NCG][j+s*N]-struct_local_ml[(NTYPE+1)*NCG+type*NCG+c][i+s*N]);
                            mean_rest[type*NCG+c] += w_inherent*diff*diff;
                        }
                    }
                    mean_den_inherent[NTYPE*NCG+c] += w_inherent;

                    double diff = (struct_local_ml[(NTYPE+1)*NCG][j+s*N]-struct_local_ml[(NTYPE+1)*NCG+NTYPE*NCG+c][i+s*N]);
                    mean_rest[NTYPE*NCG+c] += w_inherent*diff*diff;

                }
            }

            for (int type=0; type<(NTYPE+1); type++) {
                struct_local_ml[3*(NTYPE+1)*NCG+type*NCG][i+s*N] = 0.0;

            }
            for (int c=1; c<NCG; c++) {
                for (int type=0; type<(NTYPE+1); type++) {
                    struct_local_ml[3*(NTYPE+1)*NCG+type*NCG+c][i+s*N] = mean_rest[type*NCG+c]/mean_den_inherent[type*NCG+c];
                }
            }
        }
    }

    // normalize physical structural descriptors to have mean zero and unit variance
    double ** struct_local_norm; 
    struct_local_norm = dmatrix(0,N*NS-1,0,4*(NTYPE+1)*NCG-1);
    struct_mean_var = dmatrix(0,4*(NTYPE+1)*NCG-1,0,2*NTYPE-1);
    for (int k=0; k<4*(NTYPE+1)*NCG;k++) {

        double mean[NTYPE];
        double var[NTYPE];
        for (int type=0; type<NTYPE; type++) {
            mean[type] = 0.0;
            var[type] = 0.0;
        }
        //std::cout << k << std::endl;
        for (int i=0; i<N*NS;i++) {
            int type = type_data[i];
            mean[type] += struct_local_ml[k][i];
            var[type] += struct_local_ml[k][i]*struct_local_ml[k][i];
        }

        for (int type=0; type<(NTYPE); type++) {
            mean[type]/=NPerType[type]*NS;
            var[type]=sqrt(var[type]/(NPerType[type]*NS)- mean[type]*mean[type]);
            struct_mean_var[k][2*type] = mean[type];
            struct_mean_var[k][2*type+1] = var[type];
        }

        for (int i=0; i<N*NS;i++) {
            struct_local_norm[i][k] = struct_local_ml[k][i];
            //struct_local_norm[i][k] = struct_local_ml[k][i] - mean;
            //if (var > 0.00000001) struct_local_norm[i][k] /= var;
            //if (i<5 && k==NCG) std::cout << struct_local_norm[i][NCG] << std::endl;
        }
    }

    for (int type=0; type<(NTYPE); type++) {
        QString pathOrig = QString::fromStdString(folderOut);
        QString pathPred = pathOrig;
        pathPred.append("ml_struct_type");
        pathPred.append(QString("%1.csv").arg(type));
        QFile outfilePred2(pathPred);   // input file with xyz
        outfilePred2.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred2(&outfilePred2);
        //write header
        for (int c=0; c<NCG; c++) outPred2 << "DEN_CGP_TYPE0" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "DEN_CGP_TYPE1" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "DEN_CGP_TYPE2" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "DEN_CGP_ALL" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_CGP_TYPE0" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_CGP_TYPE1" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_CGP_TYPE2" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_CGP_ALL" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "PERI_CGP_TYPE0" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "PERI_CGP_TYPE1" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "PERI_CGP_TYPE2" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "PERI_CGP_ALL" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_VAR2CGP_TYPE0" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_VAR2CGP_TYPE1" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_VAR2CGP_TYPE2" << c << ",";
        for (int c=0; c<NCG; c++) outPred2 << "EPOT_VAR2CGP_ALL" << c << ",";
        for (int d=0; d<dim; d++) outPred2 << "DIM" << d << ",";
        outPred2 << "\n";

        // write body

        // length scale
        for (int k=0; k<4*(NTYPE+1)*NCG; k++) {
            outPred2 << (k % NCG)/2.0 << ",";
        }
        for (int d=0; d<dim; d++) outPred2 << "0.0" << ",";
        outPred2 << "\n";

        // mean and variance
        for (int k=0; k<4*(NTYPE+1)*NCG; k++) {
            outPred2 << struct_mean_var[k][2*type] << ",";
        }
        for (int d=0; d<dim; d++) outPred2 << "0.0" << ",";
        outPred2 << "\n";
        for (int k=0; k<4*(NTYPE+1)*NCG; k++) {
            outPred2 << struct_mean_var[k][2*type+1] << ",";
        }
        for (int d=0; d<dim; d++) outPred2 << "0.0" << ",";
        outPred2 << "\n";

        // particle data
        for (int i=0; i<N*NS; i++) {
            if (type_data[i] == type) {
                for (int k=0; k<4*(NTYPE+1)*NCG; k++) outPred2 << struct_local_norm[i][k] << ",";
                for (int d=0; d<dim; d++) outPred2 << xyz_data[i][d] << ",";
                outPred2 << "\n";
            }
        }
        outfilePred2.close();
    }
    
    free_dmatrix(struct_local_norm,0,N*NS-1,0,4*(NTYPE+1)*NCG-1);
    
}

void write_descriptors_csv_dyn(){

    // normalize physical structural descriptors to have mean zero and unit variance
    double ** dyn_mean_var; 
    dyn_mean_var = dmatrix(0,(NT+1)*NDynTotal-1,0,2*(NTYPE)-1);
    for (int k=0; k<(NT+1)*NDynTotal;k++) {

        double mean[NTYPE];
        double var[NTYPE];
        for (int type=0; type<NTYPE; type++) {
            mean[type] = 0.0;
            var[type] = 0.0;
        }
        //std::cout << k << std::endl;
        for (int i=0; i<N*NS;i++) {
            int type = type_data[i];
            if ( (k % (NT + 1) ) == NT && dyn_avg_save[i][k] > 0.000001 ) dyn_avg_save[i][k] = log10(dyn_avg_save[i][k]);
            mean[type] += dyn_avg_save[i][k];
            var[type] += dyn_avg_save[i][k]*dyn_avg_save[i][k];
        }

        for (int type=0; type<NTYPE; type++) {
            mean[type]/=NPerType[type]*NS;
            var[type]=sqrt(var[type]/(NPerType[type]*NS)- mean[type]*mean[type]);
            dyn_mean_var[k][2*type] = mean[type];
            dyn_mean_var[k][2*type+1] = var[type];
        }
    }

    for (int type=0; type<NTYPE; type++) {
        QString pathOrig = QString::fromStdString(folderOut);
        QString pathPred = pathOrig;
        pathPred.append("ml_labels_type");
        pathPred.append(QString("%1.csv").arg(type));
        QFile outfilePred(pathPred);   // input file with xyz
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred(&outfilePred);
        //write header
        for (int k=0; k<NDynTotal; k++) for (int t=1; t<=NT; t++) outPred << QString::fromStdString(DynNames[k]) << t << ",";
        outPred << "\n";

        // mean and variance
        for (int k=0; k<(NT+1)*NDynTotal; k++) {
            if ((k % (NT + 1) ) != 0) outPred << dyn_mean_var[k][2*type] << ",";
        }
        outPred<< "\n";
        for (int k=0; k<(NT+1)*NDynTotal; k++) {
            if ((k % (NT + 1) ) != 0) outPred << dyn_mean_var[k][2*type+1] << ",";
        }
        outPred << "\n";

        // write body
        for (int i=0; i<N*NS; i++) {
            if (type_data[i] == type) {
                for (int k=0; k<NDynTotal; k++) {
                    for (int t=1; t<=NT; t++) {
                        outPred << dyn_avg_save[i][k*(NT+1)+t] << ",";
                    }
                }
                outPred << "\n";

            }
        }
        outfilePred.close();
    }
    
}
