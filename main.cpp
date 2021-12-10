 
#include <QTextStream>
#include <QFile>
#include <QDir>
#include <QDebug>
#include <QString>
#include <iostream>

#include "read_lammps.h" 
#include "defs.h"
 
 int main() {
   
    // open input file
    QString path = QDir::current().path();
    path.append("/input.opt");
    QFile infile(path);   // input file with options
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);
    
    // read input file
    QString line;
    const QRegExp rxInt(QLatin1String("[^0-9]+"));
    const QRegExp rxString(QLatin1String("[^A-Z0-9]+"));
    int line_counter = 0;
    while (!infileStream.atEnd()){
        line = infileStream.readLine();
        if( line.at( 0 ) == '#' ) continue;
        switch (line_counter) {
            case 0 :
                lammpsIn = line.toStdString();
                break;
            case 1 :
                xyzOut = line.toStdString();
                break;
            case 2 : {
                 const auto&& parts = line.split(rxInt, QString::SkipEmptyParts);
                 CnfStart = parts[0].toInt();
                 CnfStep = parts[1].toInt();
                 break;
            }
            case 3 : {
                 const auto&& parts = line.split(rxInt, QString::SkipEmptyParts);
                 NS = parts[0].toInt();
                 NI = parts[1].toInt();
                 break;
            }
            case 4: {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                NDyn = parts[0].toInt();

                bb_flag = 0;
                for (int i=0 ; i<NDyn ; i++) {
                    QString obs = parts[i+1];
                    if (obs == "BB") {
                        bb_flag = 1;
                    }
                }
                break;
            }
        }
        line_counter++;
    }

    // write options
    std::cout << "Start Data Evaluation: ISOCONFIGURATIONAL ENSEMBLE" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "LAMMPS XYZ IN: " << lammpsIn << std::endl;
    std::cout << "EXTENDED XYZ OUT (DYN): " << xyzOut << std::endl;
    std::cout << "CnfStart/CnfStep: " << CnfStart << " " << CnfStep << std::endl;
    std::cout << "NS/NI: " << NS << " " << NI << std::endl;
    std::cout << "EVALUATE " << NDyn << " dynamical observables:" << std::endl;
    if (bb_flag==1) std::cout << "    Bond-Breaking" << std::endl;

 }
