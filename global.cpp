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
          if (type==0) std::cout << Ciso[type*NT+t] << std::endl;
        }
    }
    // then eval tau_alpha from isf
    for (int type=0; type<NTYPE; type++) {
      for (int t=1; t<NT; t++) {
        if (Ciso[type*NT+t] < 0.36787944117 ) {
          global_properties[type] = time_data[t]*timestep - (Ciso[type*NT+t] - 0.36787944117) * (time_data[t]*timestep-time_data[t-1]*timestep)/(Ciso[type*NT+t]-Ciso[type*NT+t-1]);
        }
      }
      if (global_properties[type] < 0.00001) { // in case simulation did not reach tau_alpha: extrapolate exponentially
        double B = -log(Ciso[type*NT+NT-1]/Ciso[type*NT+NT-2])/(time_data[NT-1]*timestep-time_data[NT-2]*timestep);
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
        global_properties[type_data[i+s*N]+NTYPE] += struct_base_local_theta_tanaka[i+s*N][0];
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
}


void calc_bonds_histograms_structure(){
    // reset histograms
    for (int j = 0; j < NCG*NStructTotal; j++) {
        for (int k = 0; k < NTYPE*NHistoStruct; k++) struct_hist[j][k] = 0.0;
        struct_ranges[j][0] = 10000;          // minimum
        struct_ranges[j][1] = -10000;         // maximum
    }

    // determine bonds and calculate hitograms
    double struct_loc[NS*N];
    for (int j = 0; j < NStructTotal; j++) {
      for (int c = 0; c < NCG; c++) {
        struct_array(j,c,struct_loc);
        for (int i=0; i<N*NS; i++) { // find minima and maxima
          if (struct_loc[i] < struct_ranges[j*NCG+c][0]) struct_ranges[j*NCG+c][0] = struct_loc[i];
          if (struct_loc[i] > struct_ranges[j*NCG+c][1]) struct_ranges[j*NCG+c][1] = struct_loc[i];
        }

        // calculate histograms
        for (int i=0; i<N*NS; i++) {
          int valint;
          if (struct_loc[i] > struct_ranges[j*NCG+c][1] - EPS) valint = NHistoStruct - 1;
          else valint = (struct_loc[i]-struct_ranges[j*NCG+c][0])/(struct_ranges[j*NCG+c][1] - struct_ranges[j*NCG+c][0])* ((double)NHistoStruct);
          struct_hist[j*NCG+c][type_data[i]*NHistoStruct + valint] += 1.0;
        }

        // normalize histograms
        for (int type=0; type<NTYPE; type++) {
          for (int k=0; k<NHistoStruct; k++) struct_hist[j*NCG+c][type*NHistoStruct+k] /= (double) NS*NPerType[type];
        }
      }
    }

}


// help functions to implement scalable code for structural descriptors
void struct_array(int index, int c, double * array){
  if (struct_base_flag >=0) {
    if (index == 0) copy_array(struct_base_local_den,c,array);
    if (index == 1) copy_array(struct_base_local_epot,c,array);
    if (index == 2) copy_array(struct_base_local_psi,2*c,array);
    if (index == 3) copy_array(struct_base_local_psi,2*c+1,array);
    if (index == 4) copy_array(struct_base_local_theta_tanaka,c,array);
    index -= 5;
  }

}
void copy_array(double ** array_in, int index, double * array_out) {
  for (int i=0; i<N*NS; i++) array_out[i]=array_in[i][index];
}