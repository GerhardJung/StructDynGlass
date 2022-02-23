
#include "eval_struct.h"
#include "defs.h"

void eval_struct(int flag){
    double Chi4[NTYPE]={}, Chi4Norm[NTYPE]={};
    std::string loc=StructNames[flag];
    loc.resize(loc.size() - 2);

    for (int type=0; type<NTYPE; type++) {
        for (int s = 0; s < NS; s++) {
            double struct_sys_avg = 0.0;
            for (int i=0; i<N;i++)
                if (type==type_data[i+s*N]) {
                    if(loc.substr( loc.length() - 2 )=="MD") {if(struct_local[NCG*(flag)][i+s*N] > 0.438) struct_sys_avg += 1.0;}
                    else if(loc.length() > 7 && loc.substr( loc.length() - 7 )=="LOG(MD)") {if(struct_local[NCG*(flag)][i+s*N] > -0.3) struct_sys_avg += 1.0;}
                    else if(loc.length() > 9 && loc.substr( loc.length() - 9 )=="LOG(FRES)") {if(struct_local[NCG*(flag)][i+s*N] > 1.0) struct_sys_avg += 1.0;}
                    else if(loc.length() > 8 && loc.substr( loc.length() - 8 )=="LOG(UTH)") {if(struct_local[NCG*(flag)][i+s*N] > -0.3) struct_sys_avg += 1.0;}
                    else {struct_sys_avg += struct_local[NCG*(flag)][i+s*N];}
                }

            Chi4Norm[type]+=struct_sys_avg;
            Chi4[type]+=struct_sys_avg*struct_sys_avg;
            //std::cout << StructNames[flag] << " " << struct_sys_avg << std::endl;
        }
        Chi4Norm[type] /= (double) NS;
        Chi4[type] /= (double) NS;
    }

    // print chi4
    QString pathPred = QString::fromStdString(folderOut);
    int t=stoi(StructNames[flag].substr( StructNames[flag].length() - 2 ));
    pathPred.append(QString("/structure_%1.dat").arg(QString::fromStdString(loc)));
    if (t==11) { // delete file
        QFile outfilePred(pathPred);   // input file with xyz
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        outfilePred.close();
    }

    QFile outfilePred2(pathPred);   // input file with xyz
    outfilePred2.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append);
    QTextStream outPred(&outfilePred2);
    outPred << time_data[t]*timestep << " ";
    for (int type=0; type<NTYPE; type++) {
            outPred << (Chi4[type] - Chi4Norm[type]*Chi4Norm[type])/((double) NPerType[type]) << " ";
    }
    outPred << "\n";
    
    outfilePred2.close();
}

void calc_bonds_histograms_structure(){
    std::cout << "CALC BONDS: STARTED" << std::endl;
    // reset histograms
    for (int j = 0; j < NCG*NStructTotal; j++) {
        for (int k = 0; k < NTYPE*NHistoStruct; k++) struct_hist[j][k] = 0.0;
        struct_ranges[j][0] = 10000;          // minimum
        struct_ranges[j][1] = -10000;         // maximum
    }

    // determine bonds and calculate hitograms
    double * struct_loc;
    for (int j = 0; j < NStructTotal; j++) {
      for (int c = 0; c < NCG; c++) {
        struct_loc = struct_local[j*NCG+c];
        //if (j==5 && c==2) std::cout << struct_ranges[j*NCG+c][0] << " " << struct_ranges[j*NCG+c][1] << " " << struct_loc[2] << std::endl;
        calc_bonds(struct_loc,struct_ranges[j*NCG+c]);
        //std::cout << StructNames[j] << " " << c << " "<< struct_ranges[j*NCG+c][0] << " " << struct_ranges[j*NCG+c][1] << std::endl;

        // calculate histograms
        for (int i=0; i<N*NS; i++) {
          int valint;
          if (struct_loc[i] > struct_ranges[j*NCG+c][1] - EPS) valint = NHistoStruct - 1;
          else valint = (struct_loc[i]-struct_ranges[j*NCG+c][0])/(struct_ranges[j*NCG+c][1] - struct_ranges[j*NCG+c][0])* ((double)NHistoStruct);
          //if (valint < 0) std::cout << struct_loc[i] << std::endl;
          struct_hist[j*NCG+c][type_data[i]*NHistoStruct + valint] += 1.0;
        }

        // normalize histograms
        for (int type=0; type<NTYPE; type++) {
          for (int k=0; k<NHistoStruct; k++) struct_hist[j*NCG+c][type*NHistoStruct+k] /= (double) NS*NPerType[type];
        }
      }
    }

      std::cout << "CALC BONDS: FINISHED " << std::endl; 
}


// help functions
void calc_bonds(double * input, double * output){
  output[0] = 10000.0;
  output[1] = -10000.0;
  for (int i=0; i<N*NS; i++) { // find minima and maxima
    if (input[i] < output[0]) output[0] = input[i];
    if (input[i] > output[1]) output[1] = input[i];
  }
}