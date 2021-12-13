
#include "defs.h"

void read_files_lammps(){

    // open first file to detect particle number, box size, dimension and timesteps
    QString path = QString::fromStdString(lammpsIn);
    path.append("/Cnf-");
    path.append(QString("%1").arg(CnfStart));
    path.append("/1.xyz");
    //std::cout << path.toStdString() << std::endl;
    QFile infile(path);   // input file with xyz
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);

    int Nloc = 0;
    QString line;
    const QRegExp rxString(QLatin1String("[^A-Z0-9]+"));
    while (!infileStream.atEnd()){
        Nloc++;
        if (!(Nloc==6) && !(Nloc==8)) line = infileStream.readLine();

        if (Nloc==4) N = line.toInt();  // particle number
        if (Nloc==6) {
            double dummy1, dummy2;
            infileStream >> dummy1 >> dummy2;
            boxL = 2.0*dummy2;
        }
        if (Nloc==8) {
            double dummy1, dummy2;
            infileStream >> dummy1 >> dummy2;
            double boxLz = 2.0*dummy2;
            if (boxLz > 1.0) dim = 3;
            else dim=2;
        }
    }
    infile.close();
    NT = Nloc/(9 + N);

    std::cout << "Simulation settings:" << std::endl;
    std::cout << "N/dim/boxL/NT: " << N << " " << dim << " " << boxL << " " << NT << std::endl;

    // allocate storage
    allocate_storage();
    
    // read all data
    path = QString::fromStdString(lammpsIn);
    path.append("/Cnf-");
    for (int i=0; i<NS; i++) {
        QString pathi = path;
        pathi.append(QString("%1/").arg(CnfStart+i*CnfStep));

        for (int j=0; j<NI; j++) {
            QString pathij = pathi;
            pathij.append(QString("%1.xyz").arg(j+1));

            QFile infile(pathij);   // input file with xyz
            infile.open(QIODevice::ReadOnly | QIODevice::Text);
            QTextStream infileStream(&infile);
            //std::cout << pathij.toStdString() << std::endl;

            for (int t=0; t<NT; t++) {
                // Skip header   
                line = infileStream.readLine();
                //if (j==0) std::cout << line.toStdString() << std::endl;
                // read time
                line = infileStream.readLine();
                time_data[t] = line.toInt();

                // skip header
                for (int l=0; l<7; l++) {
                    line = infileStream.readLine();
                }
                //if (j==0) std::cout << line.toStdString() << std::endl;
                // read data
                for (int l=0; l<N; l++) {
                    if (dim == 2) {
                        infileStream >> type_data[l+i*N] >> xyz_data[l+i*N][dim*t+dim*NT*j] >> xyz_data[l+i*N][1+dim*t+dim*NT*j];
                    } else {
                        infileStream >> type_data[l+i*N] >> xyz_data[l+i*N][dim*t+dim*NT*j] >> xyz_data[l+i*N][1+dim*t+dim*NT*j]  >> xyz_data[l+i*N][2+dim*t+dim*NT*j];
                        //if (j==0 && k==0 && l>1100) std::cout << "type " << type_data[l+i*N]<< std::endl;
                    }
                    if(i==0 && j==0 && t==0) NPerType[type_data[l+i*N]] ++;
                    
                }
                line = infileStream.readLine();
            }

            infile.close();
        }
    }

}