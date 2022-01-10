

#include "defs.h"
#include "read_write_lammps.h" 
#include "pbc.h"
#include "dyn_bb.h" 
#include "dyn_exp.h" 
#include "dyn_isf.h" 
#include "dyn_msd.h" 
#include "struct_base.h"
#include "struct_soft_modes.h"  
#include "global.h"
 
 int main() {
   
    // open input file
    QString path = QDir::current().path();
    path.append("/input.opt");
    QFile infile(path);   // input file with options
    infile.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream infileStream(&infile);
    
    // read input file
    const QRegExp rxInt(QLatin1String("[^0-9]+"));
    const QRegExp rxString(QLatin1String("[^A-Z0-9.-]+"));
    int line_counter = 0;
    int Ndyn_loc=0, Nstruct_loc=0;
    while (!infileStream.atEnd()){
        QString line = infileStream.readLine();
        if( line.at( 0 ) == '#' ) continue;
        switch (line_counter) {
            case 0 :
                lammpsIn = line.toStdString();
                break;
            case 1 :
                folderOut = line.toStdString();
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
                 NHistoStruct = parts[3].toInt();
                 break;
            }
            case 4: {
                NDyn = line.toInt();
                dyn_ranges = dmatrix(0,NDyn-1,0,1);
                DynNames = new std::string[NDyn];
                break;
            }
            case 5: {

                if (NDyn == 0) {
                    line_counter++;
                    break;
                }

                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                if (parts[0] == "BB") {
                    bb_flag = Ndyn_loc;
                    rcuti2 = parts[1].toDouble()*parts[1].toDouble();
                    rcuto2 = parts[2].toDouble()*parts[2].toDouble();
                    dyn_ranges[bb_flag][0] = parts[3].toDouble();
                    dyn_ranges[bb_flag][1] = parts[4].toDouble();
                    DynNames[bb_flag] = "BB";
                } else if (parts[0] == "EXP") {
                    exp_flag = Ndyn_loc;
                    exp_scale4i = parts[1].toDouble();
                    exp_scale4i = 1.0/(exp_scale4i*exp_scale4i*exp_scale4i*exp_scale4i);
                    dyn_ranges[exp_flag][0] = parts[2].toDouble();
                    dyn_ranges[exp_flag][1] = parts[3].toDouble();
                    DynNames[exp_flag] = "EXP";
                }  else if (parts[0] == "ISF") {
                    isf_flag = Ndyn_loc;
                    qisf = parts[1].toDouble();
                    dyn_ranges[isf_flag][0] = parts[2].toDouble();
                    dyn_ranges[isf_flag][1] = parts[3].toDouble();
                    DynNames[isf_flag] = "ISF";
                }  else if (parts[0] == "MSD") {
                    msd_flag = Ndyn_loc;
                    if (parts[2] != "DYN") {
                        dyn_ranges[msd_flag][0] = parts[2].toDouble();
                        dyn_ranges[msd_flag][1] = parts[3].toDouble();
                    } else {
                        dyn_ranges[msd_flag][0] = - 1.0; // calculate ranges dynamically
                    }
                    DynNames[msd_flag] = "MSD";
                }  

                Ndyn_loc ++;
                break;
            }
            case 6: {
                NStruct = line.toInt();
                StructNames = new std::string[NStruct+10];
                break;
            }
            case 7: {
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                if (parts[0] == "BASE") {
                    struct_base_flag = NStructTotal;
                    NStructTotal += 5;
                    NHistoGr = parts[1].toInt();
                    rcut2 = parts[2].toDouble()*parts[2].toDouble();
                    StructNames[struct_base_flag] = "DEN";
                    StructNames[struct_base_flag+1] = "EPOT";
                    StructNames[struct_base_flag+2] = "PSI5";
                    StructNames[struct_base_flag+3] = "PSI6";
                    StructNames[struct_base_flag+4] = "TT";
                }     else if (parts[0] == "SM") {
                    NStructTotal+=1;
                    struct_soft_modes_flag = NStructTotal;
                    StructNames[struct_soft_modes_flag] = "SM";
                }

                Nstruct_loc ++;
                break;
            }
        }
        line_counter++;
        if (line_counter == 6) {
            // jump back if there are still more dyn variables to come
            if (Ndyn_loc < NDyn) line_counter--;
        }
        if (line_counter == 8) {
            // jump back if there are still more dyn variables to come
            if (Nstruct_loc < NStruct) line_counter--;
        }
        
    }

    // write options
    std::cout << "Start Data Evaluation: ISOCONFIGURATIONAL ENSEMBLE" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "LAMMPS XYZ IN: " << lammpsIn << std::endl;
    std::cout << "EXTENDED XYZ OUT (DYN): " << folderOut << std::endl;
    std::cout << "CnfStart/CnfStep/timestep: " << CnfStart << " " << CnfStep << " " << timestep << std::endl;
    std::cout << "NS/NI/NHisto: " << NS << " " << NI << " " << NHisto << std::endl;
    std::cout << "EVALUATE " << NDyn << " dynamical observables:" << std::endl;
    if (bb_flag>=0) std::cout << "    Bond-Breaking" << " " << sqrt(rcuti2) << " " << sqrt(rcuto2) << std::endl;
    if (exp_flag>=0) std::cout << "    Exponential Decay" << " " << sqrt(sqrt(1.0/exp_scale4i)) << std::endl;
    if (isf_flag>=0) std::cout << "    Intermediate Scattering Function" << " " << qisf << std::endl;
    if (msd_flag>=0) std::cout << "    Mean Squared Displacement" << " " << qisf << std::endl;

    std::cout << "EVALUATE " << Nstruct_loc << " statical descriptors:" << std::endl;
    if (struct_base_flag>=0) std::cout << "    Base" << " " << NHistoGr  << std::endl;
    if (struct_soft_modes_flag>=0) std::cout << "    Soft Modes" << std::endl;

    // read files
    read_files_lammps();
    apply_pbc_global();

    // calculate statical properties
    if (struct_base_flag>=0) eval_struct_base();
    if (struct_soft_modes_flag>=0) eval_struct_soft_modes();

    // eval boundaries and structural histogramms
    calc_bonds_histograms_structure();

    // calculate isoconfigurational averages and predictabilities
    if (bb_flag>=0) eval_bb();
    if (exp_flag>=0) eval_exp();
    if (isf_flag>=0) eval_isf();
    if (msd_flag>=0) eval_msd();

    // print xyz
    print_xyz_isoconf();

    // calculate and print important global properties
    calc_print_global();
 }
