
#include "defs.h"

void read_files_lammps(){

    // open first file to detect particle number, box size, NDimension and timesteps
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
        }
    }
    infile.close();
    NT = Nloc/(9 + N);
    if (NDynTotal==0 && (struct_gnn_flag<0)) NT=1;

    std::cout << "Simulation settings:" << std::endl;
    std::cout << "N/NDim/boxL/NT: " << N << " " << NDim << " " << boxL << " " << NT << std::endl;

    // allocate storage
    allocate_storage();
    
    // read all data
    path = QString::fromStdString(lammpsIn);
    path.append("/Cnf-");
    for (int s=0; s<NS; s++) {
        QString pathi = path;
        pathi.append(QString("%1/").arg(CnfStart+s*CnfStep));

        for (int j=0; j<NI; j++) {
            QString pathij = pathi;
            QString pathij_inherent = pathi;
            pathij.append(QString("%1.xyz").arg(j+1));
            pathij_inherent.append(QString("%1_inherent.xyz").arg(j+1));

            QFile infile(pathij);   // input file with xyz
            QFile infile_inherent(pathij_inherent);   // input file with xyz inherent
            infile.open(QIODevice::ReadOnly | QIODevice::Text);
            infile_inherent.open(QIODevice::ReadOnly | QIODevice::Text);
            QTextStream infileStream(&infile);
            QTextStream infileStream_inherent(&infile_inherent);
            //std::cout << pathij.toStdString() << std::endl;

            for (int t=0; t<NT; t++) {
                // Skip header   
                line = infileStream.readLine();
                //if (j==0) std::cout << line.toStdString() << std::endl;
                // read time
                line = infileStream.readLine();
                time_data[t] = line.toInt();

                // skip header
                for (int l=0; l<7; l++) line = infileStream.readLine();
                for (int l=0; l<9; l++) line = infileStream_inherent.readLine();

                //if (j==0) std::cout << line.toStdString() << std::endl;
                // read data
                //double dmean = 0.0;
                for (int l=0; l<N; l++) {
                    double type_data_loc;
                    if (NDim == 2) {
                        infileStream >> type_data_loc >> xyz_data[l+s*N][NDim*t+NDim*NT*j] >> xyz_data[l+s*N][1+NDim*t+NDim*NT*j];
                        infileStream_inherent >> type_data_loc >> xyz_inherent_data[l+s*N][NDim*t+NDim*NT*j] >> xyz_inherent_data[l+s*N][1+NDim*t+NDim*NT*j];
                        //if (j==0 && i==0 && l==100) std::cout << "type " << type_data_loc<< " " << xyz_data[l+i*N][NDim*t+NDim*NT*j]  << " " << xyz_inherent_data[l+i*N][NDim*t+NDim*NT*j] << std::endl;
                    } else {
                        infileStream >> type_data_loc >> xyz_data[l+s*N][NDim*t+NDim*NT*j] >> xyz_data[l+s*N][1+NDim*t+NDim*NT*j]  >> xyz_data[l+s*N][2+NDim*t+NDim*NT*j];
                        infileStream_inherent >> type_data_loc >> xyz_inherent_data[l+s*N][NDim*t+NDim*NT*j] >> xyz_inherent_data[l+s*N][1+NDim*t+NDim*NT*j]  >> xyz_inherent_data[l+s*N][2+NDim*t+NDim*NT*j];
                        //if (j==0 && i==0 && l>1100) std::cout << "type " << type_data_loc<< " " << xyz_data[l+i*N][NDim*t+NDim*NT*j]  << " " << xyz_inherent_data[l+i*N][NDim*t+NDim*NT*j] << std::endl;
                    }
                    //dmean += type_data_loc;
                    if (type_cutoff[0]< 0) type_data[l+s*N] = type_data_loc - 1.0 + 0.01;
                    else {
                        int type_loc = 0;
                        while (type_data_loc > type_cutoff[type_loc] && type_loc < NTYPE -1 ) type_loc++;
                        type_data[l+s*N] = type_loc;
                        dia_data[l+s*N] = type_data_loc;
                    }
                    //std::cout << type_data_loc << " " << type_data[l+i*N] << std::endl;
                    if(s==0 && j==0 && t==0) NPerType[type_data[l+s*N]] ++;
                    
                }
                //std::cout << dmean/((double) N) << std::endl;
                line = infileStream.readLine();
                line = infileStream_inherent.readLine();
            }

            infile.close();
            infile_inherent.close();
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

    // print only first 8 files
    int NSloc = NS;
    if (NSloc > 8) NSloc = 8;
    for (int s=0; s<NSloc; s++) {
        QString pathPred = pathOrig;
        pathPred.append(QString("struct_isoconf_%1.xyz").arg(s));
        std::cout << pathPred.toStdString() << std::endl;
        QFile outfilePred(pathPred);   // input file with xyz
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred(&outfilePred);
        for (int t=1; t<=NT; t++) {
            outPred << N << "\n";
            outPred << "Properties=species:I:1:pos:R:" << NDim;
            for (int k=0; k<NStructTotal; k++) {
                if (k<struct_read_flag || k>=struct_read_flag+struct_read_Ntotal) for (int c=0; c<NCG; c++) outPred << ":" << QString::fromStdString(StructNames[k]) << c << ":R:1";
                else outPred << ":" << QString::fromStdString(StructNames[k]) << ":R:1";
            }
            for (int k=0; k<NDynTotal; k++) outPred << ":" << QString::fromStdString(DynNames[k]) << ":R:1";
            if (t<NT) outPred << " time " << time_data[t]*timestep << "\n";
            else outPred << " timescale\n";

            // calc mean
            double mean[NDynTotal] = {};
            for (int i = 0; i < N; i++) {
                for (int k=0; k<NDynTotal; k++) {
                    mean[k]+=dyn_avg_save[i+s*N][k*(NT+1)+t];
                }
            }
            if (t<NT) {
                for (int k=0; k<NDynTotal; k++) {
                    mean[k] /= (double) N;
                    for (int i = 0; i < N; i++) {
                        if (dyn_avg_save[i+s*N][k*(NT+1)+t] < mean[k]) dyn_avg_save[i+s*N][k*(NT+1)+t] = 0;
                        else dyn_avg_save[i+s*N][k*(NT+1)+t] = 1;
                    }
                }
            }

            // calc mean
            double meanS[NStructTotal] = {};
            for (int i = 0; i < N; i++) {
                for (int k=0; k<NStructTotal; k++) {
                    meanS[k]+=struct_local[k*NCG][i+s*N];
                }
            }
            for (int k=0; k<NStructTotal; k++) {
                meanS[k] /= (double) N;
                for (int i = 0; i < N; i++) {
                    if (struct_local[k*NCG][i+s*N]< meanS[k]) struct_local[k*NCG][i+s*N] = 0;
                    else struct_local[k*NCG][i+s*N] = 1;
                }
            }



            for (int i = 0; i < N; i++) {
                if (NDim == 2) {
                    outPred << type_data[i+s*N]+1 << " " << xyz_inherent_data[i+s*N][0] << " " << xyz_inherent_data[i+s*N][1] << " ";
                } else {
                    outPred << type_data[i+s*N]+1 << " " << xyz_inherent_data[i+s*N][0] << " " << xyz_inherent_data[i+s*N][1] << " " << xyz_inherent_data[i+s*N][2] << " ";
                }
                for (int k=0; k<NStructTotal; k++) {
                    if (k<struct_read_flag || k>=struct_read_flag+struct_read_Ntotal) for (int c=0; c<NCG; c++) outPred << struct_local[k*NCG+c][i+s*N] << " ";
                    else outPred << struct_local[k*NCG][i+s*N] << " ";
                }
                for (int k=0; k<NDynTotal; k++) {
                    if (t<NT) outPred << dyn_avg_save[i+s*N][k*(NT+1)+t] << " ";
                    else {
                        if (dyn_avg_save[i+s*N][k*(NT+1)+t] > 0) outPred << log10(dyn_avg_save[i+s*N][k*(NT+1)+t]) << " ";
                        else outPred << 0.0 << " ";
                    } 

                }
                outPred << "\n";
            }
        }
        outfilePred.close();
    }

    // print prediction
    if (rp_flag>=0) {
        QString pathOrig1 = QString::fromStdString(folderOut);
        QString pathPred1 = pathOrig1;
        pathPred1.append(QString("read_isoconf.dat"));
        std::cout << pathPred1.toStdString() << std::endl;
        QFile outfilePred1(pathPred1);   // input file with xyz
        outfilePred1.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred1(&outfilePred1);
        outPred1 << "TYPE";
        for (int t=1; t<NT; t++) {
            outPred1 << " LOG(FRES)" << t;
        }
        outPred1 << "\n";
        for (int i = 0; i < N*NS; i++) {
            outPred1 << type_data[i] << " ";
            for (int t=1; t<NT; t++) {
                outPred1 << dyn_avg_save[i][(rp_flag+1)*(NT+1)+t] << " ";
            }
            outPred1 << "\n";
        }
        outfilePred1.close();
    }


}