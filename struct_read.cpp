/* Implement structural descriptors as read from a specific file */

#include "struct_read.h"
#include "defs.h"
#include "pbc.h"
#include "eval_struct.h"

void eval_struct_read(){

    std::cout << "EVAL STRUCT READ " << std::endl; 

    // open first to read structural descriptor
    QString path = QString::fromStdString(struct_read_file);
    QFile infile(path);   
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);

    // read header
    QString dummy;
    infileStream >> dummy; // dummy is the type later just for double checking
    for (int k=0; k< struct_read_Ntotal; k++) {
        
         QString name;
         infileStream >> name;
         //std::cout << k << " " << name.toStdString() << std::endl;
         StructNames[struct_read_flag+k] = name.toStdString();
         //std::cout << k << " " << StructNames[struct_read_flag+k] << std::endl;
    }

    // read body
    for (int i=0; i<N*NS; i++) {
        //std::cout << i << std::endl;
        double type_test;
        infileStream >> type_test;
        int type_loc = type_test + 0.0001;
        if (type_test != type_data[i]) {
            std::cout << "Incorrect struct_read file! Line " << i << ", type(sim) " << type_data[i] << " type(ml) " << type_test << std::endl;
            //exit(0);
        }
        for (int k=0; k< struct_read_Ntotal; k++) {
           infileStream >> struct_local[NCG*(struct_read_flag+k)][i];
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