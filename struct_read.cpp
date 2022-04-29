/* Implement structural descriptors as read from a specific file */

#include "struct_read.h"
#include "defs.h"
#include "pbc.h"
#include "eval_struct.h"

void eval_struct_read(){

    std::cout << "EVAL STRUCT READ " << std::endl; 

    // open first to read structural descriptor
    QFile infiles[NTYPE];
    QTextStream infileStreams[NTYPE];

    for (int type=0; type<NTYPE; type++) {
        std::string s = std::to_string(type);
        infiles[type].setFileName(QString::fromStdString(struct_read_file + "/pred_propensity_type"+ s + ".dat"));
        infiles[type].open(QIODevice::ReadOnly | QIODevice::Text);
        infileStreams[type].setDevice(&infiles[type]);
    

        // read header
        QString dummy;
        for (int k=0; k< struct_read_Ntotal; k++) {
            
            QString name;
            infileStreams[type] >> name;
            //std::cout << k << " " << name.toStdString() << std::endl;
            StructNames[struct_read_flag+k] = name.toStdString();
            //std::cout << k << " " << StructNames[struct_read_flag+k] << std::endl;
        }
    }

    // read body
    for (int i=0; i<N*NS; i++) {
        //std::cout << i << std::endl;
        int type_loc = type_data[i];
        for (int k=0; k< struct_read_Ntotal; k++) {
           infileStreams[type_loc] >> struct_local[NCG*(struct_read_flag+k)][i];
           //if(i<5 && k==struct_read_Ntotal-1) std::cout <<struct_local[NCG*(struct_read_flag+k)][i] << std::endl;
        }
    }

    std::cout << "EVAL STRUCT READ: FINISHED1 " << std::endl; 

    // eval chi4
    for (int k=0; k< struct_read_Ntotal; k++) {
        //eval_struct(struct_read_flag+k);
    }

    std::cout << "EVAL STRUCT READ: FINISHED2 " << std::endl; 

}