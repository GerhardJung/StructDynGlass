#include "struct_gnn.h"
#include "struct_ml.h"
#include "struct_base.h"
#include "dyn_bb.h"
#include "defs.h"
#include "pbc.h"
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>


#include "ocval.h"
#include "chooseser.h"

#if defined(OC_FORCE_NAMESPACE)
using namespace OC;
#endif


void eval_struct_gnn(){

    std::cout << "EVAL GNN 1 " << std::endl; 

    // prepare structural descriptors
    // normalize physical structural descriptors to have mean zero and unit variance
    double ** struct_local_norm; 
    struct_local_norm = dmatrix(0,N*NS-1,0,5*NTYPE*NCG-1);
    for (int k=0; k<4*NTYPE*NCG;k++) {

        double mean[NTYPE];
        double var[NTYPE];
        for (int type=0; type<NTYPE; type++) {
            mean[type] = 0.0;
            var[type] = 0.0;
        }
        //std::cout << k << std::endl;
        for (int i=0; i<N*NS;i++) {
            //if (k==0 && i<1000) std::cout << struct_filion_classifiers[i][k] << std::endl;
            int type = type_data[i];
            mean[type] += struct_local_ml[k][i];
            var[type] += struct_local_ml[k][i]*struct_local_ml[k][i];
        }

        for (int type=0; type<NTYPE; type++) {
            mean[type]/=NPerType[type]*NS;
            var[type]=sqrt(var[type]/(NPerType[type]*NS)- mean[type]*mean[type]);
        }

        for (int i=0; i<N*NS;i++) {
            int type = type_data[i];
            struct_local_norm[i][k] = struct_local_ml[k][i];
            struct_local_norm[i][k] = struct_local_ml[k][i] - mean[type];
            if (var[type] > 0.00000001) struct_local_norm[i][k] /= var[type];
        }
    }

    //create folder learn and train
    QString pathOrig = QString::fromStdString(folderOut);
    pathOrig.append("/train");
    mkdir(pathOrig.toStdString().c_str(),0777);
    pathOrig = QString::fromStdString(folderOut);
    pathOrig.append("/test");
    mkdir(pathOrig.toStdString().c_str(),0777);

    for (int s=0; s<NS; s++) {

        Tab t; 

        // include box
        Arr box;
        for (int d=0; d<dim; d++) {
            box.append(boxL);
        }
        t.insertKeyAndValue("box", box);

        // include data types
        Arr types;
        for (int i=0; i<N; i++) {
            types.append(type_data[i+s*N]);
        }
        t.insertKeyAndValue("types", types);

        // include positions (t=0)
        Arr positions;
        for (int i=0; i<N; i++) {
            Arr positions_loc;
            for (int d=0; d<dim; d++) {
                positions_loc.append(xyz_data[i+s*N][d]);
            }
            positions.append(positions_loc);
        }
        t.insertKeyAndValue("positions", positions);

        // include inherent positions (t=0)
        Arr positions_inh;
        for (int i=0; i<N; i++) {
            Arr positions_loc;
            for (int d=0; d<dim; d++) {
                positions_loc.append(xyz_data[i+s*N][d]);
            }
            positions_inh.append(positions_loc);
        }
        t.insertKeyAndValue("positions_inh", positions_inh);

        // include structural descriptors
        std::string name_array[5*NTYPE] = {"DEN_CGP_TYPE0","DEN_CGP_TYPE1","DEN_CGP_TYPE2","EPOT_CGP_TYPE0","EPOT_CGP_TYPE1","EPOT_CGP_TYPE2","PERI_CGP_TYPE0",\
        "PERI_CGP_TYPE1","PERI_CGP_TYPE2","AREA_CGP_TYPE0","AREA_CGP_TYPE1","AREA_CGP_TYPE2","ASY_CGP_TYPE0","ASY_CGP_TYPE1","ASY_CGP_TYPE2"}; 
        Arr structure;
        for (int k=0; k<5*NTYPE;k++) {
            for (int cg=0; cg<NCG;cg+=2) {
                std::cout << k << " " << cg << " " << name_array[k]+"_"+std::to_string(cg) << std::endl;
                Arr structure_loc;
                for (int i=0; i<N; i++) {
                    structure_loc.append(struct_local_norm[i+s*N][k*NCG+cg]);
                }
                structure.append(structure_loc);
                t.insertKeyAndValue("structure", structure);
            }
        }
        

        // include trajectories
        //std::cout << "EVAL GNN t " << NT << std::endl; 
        /*Arr traj;
        for (int t=1; t<NT; t+=4) {
            //std::cout << "EVAL GNN t " << t << std::endl; 
            Arr traj_times;
            for (int j=0; j<NI; j++) {
                Arr traj_times_replica;
                for (int i=0; i<N; i++) {
                    Arr positions_loc;
                    for (int d=0; d<dim; d++) {
                        positions_loc.append(xyz_data[i+s*N][d+dim*t+dim*NT*j]);
                    }
                    traj_times_replica.append(positions_loc);
                }
                traj_times.append(traj_times_replica);
            }
            traj.append(traj_times);
        }
        t.insertKeyAndValue("trajectory_target_positions", traj);*/

        // include dynamical descriptors
        for (int k=0; k<NDynTotal; k++) {
            Arr dyn;
            //for (int t=1; t<=NT; t++) outPred << QString::fromStdString(DynNames[k]) << t << ",";
            for (int t=1; t<NT; t+=4) {
                //std::cout << "EVAL GNN t " << t << std::endl; 
                Arr dyn_times;
                for (int i=0; i<N; i++) {
                    dyn_times.append(dyn_avg_save[i+s*N][k*(NT+1)+t]);
                }
                dyn.append(dyn_times);
            }
            t.insertKeyAndValue(DynNames[k], dyn);
        }
        

        // write pickle
        QString pathOrig = QString::fromStdString(folderOut);
        if (s < NS/2) pathOrig.append("/train/");
        else pathOrig.append("/test/");
        pathOrig.append(QString("aggregated_%1.pickle").arg(s));
        DumpValToFile(t, pathOrig.toStdString(), SERIALIZE_P2 );

    }

}