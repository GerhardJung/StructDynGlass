

#include "defs.h"
#include "read_write_lammps.h" 
#include "pbc.h"
#include "dyn_bb.h" 
 
 int main() {
   
    // open input file
    QString path = QDir::current().path();
    path.append("/input.opt");
    QFile infile(path);   // input file with options
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);
    
    // read input file
    const QRegExp rxInt(QLatin1String("[^0-9]+"));
    const QRegExp rxString(QLatin1String("[^A-Z0-9.]+"));
    int line_counter = 0;
    int Ndyn_loc=0;
    while (!infileStream.atEnd()){
        QString line = infileStream.readLine();
        if( line.at( 0 ) == '#' ) continue;
        switch (line_counter) {
            case 0 :
                lammpsIn = line.toStdString();
                break;
            case 1 :
                xyzOut = line.toStdString();
                break;
            case 2 : {
                 const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                 CnfStart = parts[0].toInt();
                 CnfStep = parts[1].toInt();
                 timestep = parts[2].toDouble();
                 break;
            }
            case 3 : {
                 const auto&& parts = line.split(rxInt, QString::SkipEmptyParts);
                 NS = parts[0].toInt();
                 NI = parts[1].toInt();
                 NHisto = parts[2].toInt();
                 break;
            }
            case 4: {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                NDyn = line.toInt();
                break;
            }
            case 5: {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                if (parts[0] == "BB") {
                    bb_flag = 1;
                    rcuti2 = parts[1].toDouble()*parts[1].toDouble();
                    rcuto2 = parts[2].toDouble()*parts[2].toDouble();
                }        


                Ndyn_loc ++;
                break;
            }
        }
        line_counter++;
        if (line_counter == 6) {
            // jump back if there are still more dyn variables to come
            if (Ndyn_loc < NDyn) line_counter--;
        }
    }

    // write options
    std::cout << "Start Data Evaluation: ISOCONFIGURATIONAL ENSEMBLE" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "LAMMPS XYZ IN: " << lammpsIn << std::endl;
    std::cout << "EXTENDED XYZ OUT (DYN): " << xyzOut << std::endl;
    std::cout << "CnfStart/CnfStep/timestep: " << CnfStart << " " << CnfStep << " " << timestep << std::endl;
    std::cout << "NS/NI/NHisto: " << NS << " " << NI << " " << NHisto << std::endl;
    std::cout << "EVALUATE " << NDyn << " dynamical observables:" << std::endl;
    if (bb_flag==1) std::cout << "    Bond-Breaking" << " " << sqrt(rcuti2) << " " << sqrt(rcuto2) << std::endl;


    // read files
    read_files_lammps();
    //apply_pbc_global();

    // calculate isoconfigurational averages and predictabilities
    if (bb_flag) eval_bb();

 }
