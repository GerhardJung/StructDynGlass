
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
        if (!(Nloc==7) && !(Nloc==9)) line = infileStream.readLine();

        if (Nloc==4) N = line.toInt();  // particle number
        if (Nloc==7) {
            double dummy1, dummy2;
            infileStream >> dummy1 >> dummy2;
            boxL = 2.0*dummy2;
        }
        if (Nloc==9) {
            double dummy1, dummy2;
            infileStream >> dummy1 >> dummy2;
            double boxLz = 2.0*dummy2;
            //std::cout << boxL << " " << boxLz << std::endl;
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
                    type_data[l+i*N] --;
                    if(i==0 && j==0 && t==0) NPerType[type_data[l+i*N]] ++;
                    
                }
                line = infileStream.readLine();
            }

            infile.close();
        }
    }
    std::cout << "NPerType ";
    for (int type=0;type< NTYPE; type++) {
        std::cout << type << ": " <<  NPerType[type] << " ";
    }
    std::cout << std::endl;

    // calculate q value for ISF
    double qmin = 2*M_PI/boxL;
    double Nqmin = qisf / qmin;
    int Nqmin_int = Nqmin + 0.5;
    double qisf = Nqmin_int*qmin;
    std::cout << "UPDATED: qisf to fit into box. New value " << qisf << std::endl;

    /*for (int l=0; l<N; l++) {
        std::cout << type_data[l+0*N] << std::endl;
    }*/

}


// print extended xyz file to visualized isoconfigurational average
void print_xyz_isoconf(){
    QString pathOrig = QString::fromStdString(folderOut);

    // print predictabilities
    for (int s=0; s<NS; s++) {
        QString pathPred = pathOrig;
        pathPred.append(QString("/struct_isoconf_%1.xyz").arg(s));
        std::cout << pathPred.toStdString() << std::endl;
        QFile outfilePred(pathPred);   // input file with xyz
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred(&outfilePred);
        for (int t=1; t<NT; t++) {
            outPred << N << "\n";
            outPred << "Properties=species:I:1:pos:R:" << dim << ":";
            if (bb_flag) outPred << "BB:R:1";
            if (struct_base_flag) outPred << "Psi:R:2:TT:R:1";
            outPred << " time " << time_data[t]*timestep << "\n";
            for (int i = 0; i < N; i++) {
                if (dim == 2) {
                    outPred << type_data[i+s*N]+1 << " " << xyz_data[i+s*N][0] << " " << xyz_data[i+s*N][1] << " ";
                } else {
                    outPred << type_data[i+s*N]+1 << " " << xyz_data[i+s*N][0] << " " << xyz_data[i+s*N][1] << " " << xyz_data[i+s*N][2] << " ";
                }
                if (bb_flag) outPred << dyn_bb_avg[i+s*N][t] << " ";
                if (struct_base_flag) outPred << struct_base_local_psi[i+s*N][0] << " " << struct_base_local_psi[i+s*N][1] << " " << struct_base_local_theta_tanaka[i+s*N];
                outPred << "\n";
            }
        }
        outfilePred.close();
    }

}