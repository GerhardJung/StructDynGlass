#include "global.h"
#include "defs.h"


void calc_print_global(){
  // calc tau_alpha
  if (isf_flag>=0) {
    // first calc average isf
    double Ciso[NTYPE*NT];
    for (int type=0; type<NTYPE*NT; type++) {
        Ciso[type] = 0.0;
    }
    for (int s=0; s<NS;s++) { // loop over structures
      for (int i=0; i<N;i++) { // loop over particles
        for (int t=0; t<NT; t++) {
          Ciso[type_data[i+s*N]*NT+t] += dyn_avg_save[i+s*N][t+isf_flag*NT];
        }
      }
    }
    for (int type=0; type<NTYPE; type++) {
        for (int t=0; t<NT; t++) {
          Ciso[type*NT+t] /= (double) NS*NPerType[type];
          //if (type==0) std::cout << t << " " << Ciso[type*NT+t] << std::endl;
        }
    }
    // then eval tau_alpha from isf
    for (int type=0; type<NTYPE; type++) {
      for (int t=1; t<NT; t++) {
        if (Ciso[type*NT+t] < 0.36787944117 ) {
          global_properties[type] = time_data[t]*timestep - (Ciso[type*NT+t] - 0.36787944117) * (time_data[t]*timestep-time_data[t-1]*timestep)/(Ciso[type*NT+t]-Ciso[type*NT+t-1]);
          break;
          //std::cout << global_properties[type] << std::endl;
        }
      }
      if (global_properties[type] < 0.00001) { // in case simulation did not reach tau_alpha: extrapolate exponentially
        double B = -log(Ciso[type*NT+NT-1]/Ciso[type*NT+NT-3])/(time_data[NT-1]*timestep-time_data[NT-3]*timestep);
        double A = Ciso[type*NT+NT-1]/exp(-B*time_data[NT-1]*timestep);
        //std::cout << A << " " << B << std::endl;
        global_properties[type] = -log(0.36787944117/A)/B;
      }
    }
  }

  // calc global theta tanaka
  if (struct_base_flag>=0) {
    for (int s=0; s<NS;s++) { // loop over structures
      for (int i=0; i<N;i++) { // loop over particles
        global_properties[type_data[i+s*N]+NTYPE] += struct_local[NCG*(struct_base_flag+4)][i+s*N];
      }
    }
  }
  for (int type=0; type<NTYPE; type++) {
     global_properties[type+NTYPE] /= (double) NS*NPerType[type];
  }


  // print global properties
  QString pathGlobal = QString::fromStdString(folderOut);
  pathGlobal.append("/global_properties.dat");
  QFile outfileGlobal(pathGlobal);   // input file with xyz
  outfileGlobal.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream outGlobal(&outfileGlobal);
  outGlobal << "parameter ";
  for (int type=0; type<NTYPE; type++) {
    outGlobal << type+1 << " ";
  }
  outGlobal << "\n";
  for (int k=0; k<2; k++) {
      if (k==0) outGlobal << "tau_alpha ";
      if (k==1) outGlobal << "theta_tanaka ";
      for (int type=0; type<NTYPE; type++) {
        outGlobal << global_properties[type+k*NTYPE] << " ";
      }
      outGlobal << "\n";
  }
  outfileGlobal.close();


  // calculate soft mode plots
  if (struct_soft_modes_flag >= 0) {
     // find wmax
     double wmax = 0.0;
     for (int s=0; s<NS;s++) {
        if( hessian_evalues[s*N*dim+N*dim-dim-1] > wmax ) wmax = hessian_evalues[s*N*dim+N*dim-dim-1];
     }

     // reset histogram
     for (int l=0; l<NHistoSM; l++) {
       sm_histograms[l][0] = sm_histograms[l][1] = sm_histograms[l][2] = 0.0;
     }

     // fill histograms
     int Nquasi = 0;
     for (int s=0; s<NS;s++) {
        for (int k=0; k<N*dim-dim;k++) {
            double Pw = 0.0;
            for (int i=0; i<N;i++) {
                double resloc=0.0;
                for (int di=0; di<dim;di++) resloc+=hessian_evectors[s*N*dim+k][i*dim+di]*hessian_evectors[s*N*dim+k][i*dim+di];
                Pw += resloc*resloc;
            }
            
            int wloc = hessian_evalues[s*N*dim+k]/wmax*NHistoSM;
            if (wloc < NHistoSM && wloc >=0) {
              sm_histograms[wloc][0] += 1.0;
              if (1.0/(N*Pw) < Pcut) {
                sm_histograms[wloc][1] += 1.0;
                Nquasi ++;
              }
              sm_histograms[wloc][2] += 1.0/(N*Pw);
            }
        }
    }

    // print histograms
    pathGlobal = QString::fromStdString(folderOut);
    pathGlobal.append("/sm_histograms.dat");
    QFile outfileGlobal2(pathGlobal);   // input file with xyz
    outfileGlobal2.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outGlobal2(&outfileGlobal2);
    outGlobal2 << "time D(w) D_loc(w) P(w)\n";
    for (int l=0; l<NHistoSM; l++) {
        if (sm_histograms[l][0] > 0) {
          outGlobal2 << (l+0.5)*wmax/((double)NHistoSM) << " " << sm_histograms[l][0]/((double) (dim*N-dim)*NS)<< " " << sm_histograms[l][1]/((double) Nquasi)<< " " << sm_histograms[l][2]/sm_histograms[l][0] ;
          outGlobal2 << "\n";
        }
    }
    outfileGlobal2.close();

    pathGlobal = QString::fromStdString(folderOut);
    pathGlobal.append("/sm_participation.dat");
    QFile outfileGlobal3(pathGlobal);   // input file with xyz
    outfileGlobal3.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outGlobal3(&outfileGlobal3);
    outGlobal3 << "w P(w)\n";
     for (int s=0; s<NS;s++) {
        for (int k=0; k<N*dim-dim;k++) {
            double Pw = 0.0;
            for (int i=0; i<N;i++) {
                double resloc=0.0;
                for (int di=0; di<dim;di++) resloc+=hessian_evectors[s*N*dim+k][i*dim+di]*hessian_evectors[s*N*dim+k][i*dim+di];
                Pw += resloc*resloc;
            }
            if (hessian_evalues[s*N*dim+k] > 0.0) outGlobal3 << hessian_evalues[s*N*dim+k]<< " " << 1.0/(N*Pw) << "\n";
        }
    }
    outfileGlobal.close();

  }
}


void calc_bonds_histograms_structure(){
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
        //if (j==5 && c==2) std::cout << struct_ranges[j*NCG+c][0] << " " << struct_ranges[j*NCG+c][1] << std::endl;

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