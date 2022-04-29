#include "global.h"
#include "defs.h"

int Nglobal=0;

void calc_print_global(){
  // calc tau_alpha
  double tau_alpha_isf;
  if (isf_flag>=0) {
     calc_tau_alpha(isf_flag,Nglobal);
     Nglobal++;
  }
  if (exp_flag>=0) {
     calc_tau_alpha(exp_flag,Nglobal);
     Nglobal++;
  }
  if (bb_flag>=0) {
     calc_tau_alpha(bb_flag,Nglobal);
     Nglobal++;
  }
  if (rp_flag>=0) {
     calc_tau_alpha(rp_flag+2,Nglobal);
     Nglobal++;
  }

  // calc global theta tanaka
  if (struct_base_flag>=0) {
    for (int s=0; s<NS;s++) { // loop over structures
      for (int i=0; i<N;i++) { // loop over particles
        global_properties[type_data[i+s*N]+Nglobal*NTYPE] += struct_local[NCG*(struct_base_flag+4)][i+s*N];
      }
    }
  }
  for (int type=0; type<NTYPE; type++) {
     global_properties[type+Nglobal*NTYPE] /= (double) NS*NPerType[type];
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
  int k=0;
  if (isf_flag>=0) { 
    outGlobal << "tau_alpha_isf ";
    for (int type=0; type<NTYPE; type++) {
      outGlobal << global_properties[type+k*NTYPE] << " ";
    }
    outGlobal << "\n";
    k++;
  }
  if (exp_flag>=0) { 
      outGlobal << "tau_alpha_exp ";
    for (int type=0; type<NTYPE; type++) {
      outGlobal << global_properties[type+k*NTYPE] << " ";
    }
    outGlobal << "\n";
    k++;
  }
  if (bb_flag>=0) { 
      outGlobal << "tau_alpha_bb ";
    for (int type=0; type<NTYPE; type++) {
      outGlobal << global_properties[type+k*NTYPE] << " ";
    }
    outGlobal << "\n";
    k++;
  }
  if (rp_flag>=0) { 
      outGlobal << "tau_alpha_rp ";
    for (int type=0; type<NTYPE; type++) {
      outGlobal << global_properties[type+k*NTYPE] << " ";
    }
    outGlobal << "\n";
    k++;
  }
  if (struct_base_flag>=0) { 
    outGlobal << "theta_tanaka ";
    for (int type=0; type<NTYPE; type++) {
      outGlobal << global_properties[type+k*NTYPE] << " ";
    }
    outGlobal << "\n";
    k++;
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
            double Pw = 1.0/(N*participation_ratio[s*N*dim+k]);
            
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

void calc_tau_alpha(int flag, int count){
    double Ciso[NTYPE*NT];
    for (int type=0; type<NTYPE*NT; type++) {
        Ciso[type] = 0.0;
    }
    for (int s=0; s<NS;s++) { // loop over structures
      for (int i=0; i<N;i++) { // loop over particles
        for (int t=0; t<NT; t++) {
          Ciso[type_data[i+s*N]*NT+t] += dyn_avg_save[i+s*N][t+flag*(NT+1)];
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
          global_properties[count*NTYPE+type] = time_data[t]*timestep - (Ciso[type*NT+t] - 0.36787944117) * (time_data[t]*timestep-time_data[t-1]*timestep)/(Ciso[type*NT+t]-Ciso[type*NT+t-1]);
          break;
          //std::cout << global_properties[type] << std::endl;
        }
      }
      if (global_properties[count*NTYPE+type] < 0.00001) { // in case simulation did not reach tau_alpha: extrapolate exponentially
        double B = -log(Ciso[type*NT+NT-1]/Ciso[type*NT+NT-3])/(time_data[NT-1]*timestep-time_data[NT-3]*timestep);
        double A = Ciso[type*NT+NT-1]/exp(-B*time_data[NT-1]*timestep);
        //std::cout << " " << flag << " " << A << " " << B << std::endl;
        global_properties[count*NTYPE+type] = -log(0.36787944117/A)/B;
      }
    }
}
