

#include "defs.h"
#include "read_write_lammps.h" 
#include "pbc.h"
#include "dyn_bb.h" 
#include "dyn_exp.h" 
#include "dyn_isf.h" 
#include "dyn_msd.h" 
#include "dyn_rearrange_patinet.h" 
#include "struct_base.h"
#include "struct_soft_modes.h"  
#include "struct_filion.h"  
#include "struct_voronoi.h"  
#include "struct_read.h"  
#include "struct_epot.h" 
#include "struct_ml.h"
#include "struct_gnn.h"
#include "eval_struct.h"
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
    const QRegExp rxString(QLatin1String("[^A-Za-z0-9./_-]+"));
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

            case 4 : {
                 const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                 int NDim_loc = parts[0].toInt();
                 if (NDim_loc != NDim) {
                     printf("Program compiled with NDimension %d but input script requests NDimension %d!\n",NDim,NDim_loc);
                     exit(0);
                 }
                 model = parts[1].toStdString();
                 NTYPE = parts[2].toInt();
                 type_cutoff= dvector(0,NTYPE-2); 
                 if (parts[3] == "DISCRETE") {
                     type_cutoff[0] = -1.0;
                 } else {
                     for (int type=0;type<NTYPE-1; type++) {
                         type_cutoff[type] = parts[3+type].toDouble();
                     }
                 }
                 break;
            }


            case 5: {
                NDyn = line.toInt();
                dyn_ranges = dmatrix(0,NDyn+40-1,0,2);
                DynNames = new std::string[NDyn+40];
                if (NDyn == 0) {
                    line_counter++;
                }
                break;
            }
            case 6: {

                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                if (parts[0] == "BB") {
                    bb_flag = NDynTotal;
                    NDynTotal+=1;
#ifdef USE_RELATIVE
                    NDynTotal+=N_NEIGH_MAX;
#endif
                    rcuti2 = parts[1].toDouble()*parts[1].toDouble();
                    rcuto2 = parts[2].toDouble()*parts[2].toDouble();
                    dyn_ranges[bb_flag][0] = parts[3].toDouble();
                    dyn_ranges[bb_flag][1] = parts[4].toDouble();
                    DynNames[bb_flag] = "BB";
#ifdef USE_RELATIVE
                    for (int i=0; i<N_NEIGH_MAX; i++)  {
                        DynNames[bb_flag+i+1] = "DRREL" + std::to_string(i) + "X";
                        dyn_ranges[bb_flag+i+1][0] = 0.0;
                        dyn_ranges[bb_flag+i+1][1] = 1.0;
                    }
#endif
                } else if (parts[0] == "EXP") {
                    exp_flag = NDynTotal;
                    NDynTotal+=1;
                    exp_scale4i = parts[1].toDouble();
                    exp_scale4i = 1.0/(exp_scale4i*exp_scale4i*exp_scale4i*exp_scale4i);
                    dyn_ranges[exp_flag][0] = parts[2].toDouble();
                    dyn_ranges[exp_flag][1] = parts[3].toDouble();
                    DynNames[exp_flag] = "EXP";
                }  else if (parts[0] == "ISF") {
                    isf_flag = NDynTotal;
                    NDynTotal+=1;
                    qisf = parts[1].toDouble();
                    dyn_ranges[isf_flag][0] = parts[2].toDouble();
                    dyn_ranges[isf_flag][1] = parts[3].toDouble();
                    DynNames[isf_flag] = "ISF";
                }  else if (parts[0] == "MSD") {
                    msd_flag = NDynTotal;
                    NDynTotal+=2;
                    if (parts[1] != "DYN") {
                        dyn_ranges[msd_flag][0] = parts[1].toDouble();
                        dyn_ranges[msd_flag][1] = parts[2].toDouble();
                    } else {
                        dyn_ranges[msd_flag][0] = - 10001.0; // calculate ranges dynamically
                    }
                    if (parts[1] != "DYN") {
                        dyn_ranges[msd_flag+1][0] = parts[3].toDouble();
                        dyn_ranges[msd_flag+1][1] = parts[4].toDouble();
                    } else {
                        dyn_ranges[msd_flag+1][0] = - 10001.0; // calculate ranges dynamically
                    }
                    if (parts.size() > 5) overlap_cut = parts[5].toDouble();
                    DynNames[msd_flag] = "MD";
                    DynNames[msd_flag+1] = "LOG(MD)";
                }  else if (parts[0] == "RP") {
                    inherent = 1;
                    rp_flag = NDynTotal;
                    NDynTotal+=3;
                    dyn_rearrange_mode = parts[1].toInt();
                    dyn_rearrange_threshold = parts[2].toDouble();
                    if (parts[2] != "DYN") {
                        dyn_ranges[rp_flag][0] = parts[3].toDouble();
                        dyn_ranges[rp_flag][1] = parts[4].toDouble();
                    } else {
                        dyn_ranges[rp_flag][0] = - 10001.0; // calculate ranges dynamically
                    }
                    if (parts[4] != "DYN") {
                        dyn_ranges[rp_flag+1][0] = parts[5].toDouble();
                        dyn_ranges[rp_flag+1][1] = parts[6].toDouble();
                    } else {
                        dyn_ranges[rp_flag+1][0] = - 10001.0; // calculate ranges dynamically
                    }
                    if (parts[6] != "DYN") {
                        dyn_ranges[rp_flag+2][0] = parts[7].toDouble();
                        dyn_ranges[rp_flag+2][1] = parts[8].toDouble();
                    } else {
                        dyn_ranges[rp_flag+2][0] = - 10001.0; // calculate ranges dynamically
                    }
                    DynNames[rp_flag] = "LOG(UTH)";
                    DynNames[rp_flag+1] = "LOG(FRES)";
                    DynNames[rp_flag+2] = "RP";
                }  

                Ndyn_loc ++;
                break;
            }
            case 7: {
                NStruct = line.toInt();
                StructNames = new std::string[NStruct+550];
                break;
            }
            case 8: {
                if (NStruct == 0) continue;
                const auto&& parts = line.split(rxString, QString::SkipEmptyParts);
                if (parts[0] == "BASE") {

                    struct_base_flag = NStructTotal;
                    NStructTotal += 6+2*(lmax-lmin);
                    NHistoGr = parts[1].toInt();
                    rcut2 = parts[2].toDouble()*parts[2].toDouble();
                    StructNames[struct_base_flag] = "DEN";
                    StructNames[struct_base_flag+2] = "EPOT";
                    StructNames[struct_base_flag+4] = "TT";
                    for (int l=lmin; l<lmax; l++) {
                        StructNames[struct_base_flag+6 + 2*l - 2*lmin] = "PSI" + std::to_string(l);
                    }
                    for (int k=0; k<3 + lmax - lmin; k++) {
                        StructNames[struct_base_flag+2*k+1] = StructNames[struct_base_flag+2*k] + "_INH";
                    }
                }     else if (parts[0] == "SM") {
                    struct_soft_modes_flag = NStructTotal;
                    NStructTotal+=2;
                    NHistoSM = parts[1].toInt();
                    Pcut = parts[2].toDouble();
                    if (parts[3] == "READ") modeSM = 1;
                    else modeSM = 0;
                    StructNames[struct_soft_modes_flag] = "SM";
                    StructNames[struct_soft_modes_flag+1] = "VIB";
                } else if (parts[0] == "MLFILION") {
                    struct_filion_flag = NStructTotal;
                    struct_filion_file = parts[1].toStdString();
                } else if (parts[0] == "MLJUNG") {
                    struct_ml_flag = NStructTotal;
                    if (parts.size() > 1) if (parts[1] == "CALCDR") calc_dr = 1;
                } else if (parts[0] == "MLGNN") {
                    struct_gnn_flag = NStructTotal;
                } else if (parts[0] == "EPOT") {
                    struct_epot_flag = NStructTotal;
                } else if (parts[0] == "VORONOI") {
                    struct_voronoi_flag = NStructTotal;
                    NStructTotal+=3;
                } else if (parts[0] == "READ") {
                    struct_read_flag = NStructTotal;
                    struct_read_Ntotal = parts[1].toInt();
                    NStructTotal += struct_read_Ntotal;
                    struct_read_file = parts[2].toStdString();
                }

                Nstruct_loc ++;
                break;
            }
        }
        line_counter++;
        //printf("%d\n",line_counter);
        if (line_counter == 7) {
            // jump back if there are still more dyn variables to come
            if (Ndyn_loc < NDyn) line_counter--;
        }
        if (line_counter == 9) {
            // jump back if there are still more dyn variables to come
            if (Nstruct_loc < NStruct) line_counter--;
        }
        
    }
    if (inherent==0 && NStruct==0) noinherent = 1;

    // write options
    std::cout << "Start Data Evaluation: ISOCONFIGURATIONAL ENSEMBLE" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "LAMMPS XYZ IN: " << lammpsIn << std::endl;
    std::cout << "EXTENDED XYZ OUT (DYN): " << folderOut << std::endl;
    std::cout << "CnfStart/CnfStep/timestep: " << CnfStart << " " << CnfStep << " " << timestep << std::endl;
    std::cout << "NS/NI/NHisto/NHistoStruct: " << NS << " " << NI << " " << NHisto << " " << NHistoStruct << std::endl;
    std::cout << "NDimension/Model/NTypes: " << NDim << " " << model << " " << NTYPE  << std::endl;
    std::cout << "EVALUATE " << NDyn << " dynamical observables:" << std::endl;
    if (bb_flag>=0) std::cout << "    " << bb_flag << ": Bond-Breaking " << sqrt(rcuti2) << " " << sqrt(rcuto2) << std::endl;
    if (exp_flag>=0) std::cout << "    " << exp_flag << ": Exponential Decay" << " " << sqrt(sqrt(1.0/exp_scale4i)) << std::endl;
    if (isf_flag>=0) std::cout << "    " << isf_flag << ": Intermediate Scattering Function" << " " << qisf << std::endl;
    if (rp_flag>=0) std::cout << "    " << rp_flag << ": Patinet Rearrangement Detection" << std::endl;
    if (msd_flag>=0) std::cout << "    " << msd_flag << ": Mean Displacement"  << std::endl;

    std::cout << "EVALUATE " << Nstruct_loc << " statical descriptors:" << std::endl;
    if (struct_base_flag>=0) std::cout << "    " << struct_base_flag << ": Base" << " " << NHistoGr  << std::endl;
    if (struct_soft_modes_flag>=0) std::cout << "    " << struct_soft_modes_flag << ": Soft Modes" << std::endl;
    if (struct_filion_flag>=0) std::cout << "    " << struct_filion_flag << ": ML Filion" << std::endl;
    if (struct_ml_flag>=0) std::cout << "    " << struct_ml_flag << ": ML Jung" << std::endl;
    if (struct_gnn_flag>=0) std::cout << "    " << struct_gnn_flag << ": ML GNN" << std::endl;
    if (struct_read_flag>=0) std::cout << "    " << struct_read_flag << ": Read Structural Descriptors " << struct_read_file << std::endl;
    if (struct_voronoi_flag>=0) std::cout << "    " << struct_voronoi_flag << ": CALC Voronoi " << std::endl;

    // read files
    read_files_lammps();
    apply_pbc_global();

    // calculate statical properties
    if (struct_base_flag>=0) eval_struct_base();
    if (struct_soft_modes_flag>=0) eval_struct_soft_modes();
    if (struct_filion_flag>=0) eval_struct_filion();
    if (struct_ml_flag>=0) eval_struct_ml();
    if (struct_read_flag>=0) eval_struct_read();
    if (struct_voronoi_flag>=0) eval_struct_voronoi();
     if (struct_epot_flag>=0) eval_struct_epot();

    // print out extended physical descriptors
    if (struct_base_flag>=0) write_descriptors_csv_phys();

    // eval boundaries and structural histogramms
    calc_bonds_histograms_structure();

    // calculate isoconfigurational averages and predictabilities
    if (bb_flag>=0) eval_bb();
    if (exp_flag>=0) eval_exp();
    if (isf_flag>=0) eval_isf();
    if (rp_flag>=0) eval_rp();
    if (msd_flag>=0) eval_msd();

    // print learning batches for machine learning
    if (NDynTotal>0) write_descriptors_csv_dyn();

    // write gnn
    if (struct_gnn_flag>=0) eval_struct_gnn();

    // print xyz
    print_xyz_isoconf();

    // calculate and print important global properties
    calc_print_global();
 }
