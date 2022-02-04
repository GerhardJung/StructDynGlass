/* Implement structural descriptors as read from a specific file */

#include "struct_read.h"
#include "defs.h"
#include "pbc.h"

void eval_struct_read(){
    // open first to read structural descriptor
    QString path = QString::fromStdString(lammpsIn);
    path.append(QString::fromStdString(struct_read_file));
    QFile infile(path);   
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);

    // read header
    QString dummy;
    infileStream >> dummy;
    for (int k=0; k< struct_read_Ntotal; k++) {
         QString name;
         infileStream >> name;
         StructNames[struct_read_flag+k] = name.toStdString();
    }

    // read body
    for (int i=0; i<N*NS; i++) {
        int type_test;
        infileStream >> type_test;
        if (type_test != type_data[i]) {
            std::cout << "Incorrect struct_read file!" << std::endl;
            exit(0);
        }
        for (int k=0; k< struct_read_Ntotal; k++) {
           infileStream >> struct_local[struct_read_flag+k][i];
        }
    }

}