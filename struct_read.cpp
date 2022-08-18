/* Implement structural descriptors as read from a specific file */

#include "struct_read.h"
#include "defs.h"
#include "pbc.h"
#include "eval_struct.h"

void eval_struct_read(){

    std::cout << "EVAL STRUCT READ " << std::endl; 

    // open first to read structural descriptor
    QFile infiles;
    QTextStream infileStreams;

    //for (int type=0; type<NTYPE; type++) {
    //    std::string s = std::to_string(type);
    infiles.setFileName(QString::fromStdString(struct_read_file + "/pred_propensity.dat"));
    infiles.open(QIODevice::ReadOnly | QIODevice::Text);
    infileStreams.setDevice(&infiles);
    
    // read header
    QString name;
    // read first column (ID)
    infileStreams >> name;
    for (int k=0; k< struct_read_Ntotal; k++) {
        
        infileStreams >> name;
        //std::cout << k << " " << name.toStdString() << std::endl;
        StructNames[struct_read_flag+k] = name.toStdString();
        //std::cout << k << " " << StructNames[struct_read_flag+k] << std::endl;
    }
    //}

    // read body
    for (int i=0; i<N*NS; i++) {
        int j;
        infileStreams >> j;
        if (j !=i) std::cout << "WRONG INDEX " << i << " " << j << std::endl;
        for (int k=0; k< struct_read_Ntotal; k++) {
           infileStreams >> struct_local[NCG*(struct_read_flag+k)][i];
           //if(i<5 && k==struct_read_Ntotal-1) std::cout <<struct_local[NCG*(struct_read_flag+k)][i] << std::endl;
        }
    }

    std::cout << "EVAL STRUCT READ: FINISHED1 " << std::endl; 

    // eval chi4, S4 and G4
    for (int k=0; k< struct_read_Ntotal; k++) {
        int first = 0;
        if (k==0) first = 1;
        eval_struct(struct_local[NCG*(struct_read_flag+k)],StructNames[struct_read_flag+k],first);
    }

    std::cout << "EVAL STRUCT READ: FINISHED2 " << std::endl; 

}