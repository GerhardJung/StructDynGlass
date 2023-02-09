 
#include "eval_isoconf.h"
#include "global.h"
#include "pbc.h"
#include <bits/stdc++.h> 

using namespace std;

void reset_dyn(int t){
    for (int i = 0; i < NS*N; i++) {
        for (int j = 0; j < NHisto; j++) {
            dyn_hist_data[i][j] = 0.0;
        }
        dyn_avg[i] = 0.0;
        dyn_avg2[i] = 0.0;
    }
    for (int i = 0; i < NS*NI*NTYPE; i++) dyn_sys_avg[i] = 0.0;
}

void add_histogram_avg(int s, int i, int j,double* hist_cut, double val){
    int valint;
    double hist_lower = hist_cut[0];
    double hist_upper = hist_cut[1];
    if (val > hist_upper - EPS && val < hist_upper + EPS) valint = NHisto - 1;
    else valint = (val-hist_lower)/(hist_upper - hist_lower)* ((double)NHisto);
    //if ( valint < 0 || valint >= NHisto) std::cout << valint << " " << val << " " <<  (val-hist_lower)/(hist_upper - hist_lower) << " " <<  hist_upper << " " <<  hist_lower << std::endl;
    //std::cout << val << " " << valint << std::endl;
    if (valint >=0 && valint < NHisto) dyn_hist_data[s*N+i][valint]+=1.0;

    dyn_avg[s*N+i] += val;
    dyn_avg2[s*N+i] += val*val; 
    //if (val < hist_cut[2]) dyn_sys_avg[s*NTYPE*NI+type_data[s*N+i]*NI+j] += 1.0;
}

void norm_histogram(){
    double norm;
    for (int i = 0; i < NS*N; i++) {
        norm = 0.0;
        for (int j = 0; j < NHisto; j++) {
            norm += dyn_hist_data[i][j];
        }
        for (int j = 0; j < NHisto; j++) {
            dyn_hist_data[i][j]/=norm;
        }
    }
}

void eval_isoconf(int t, int flag){

    norm_histogram();

    // print per particle histo
    // print histograms
    //if (t==13) {
    /*    QString pathHisto =  QString::fromStdString(folderOut);
        pathHisto.append(QString("/isoconf_check_%1_t%2.dat").arg(QString::fromStdString(dyn)).arg(t));
        QFile outfileHisto(pathHisto);   // input file with xyz
        outfileHisto.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outHisto(&outfileHisto);
        for (int i = 0; i < NS*N; i++) {
            outHisto << i << " ";
            for (int j = 0; j < NHisto; j++) {  
                outHisto << dyn_hist_data[i][j]  << " ";
            }
            outHisto << "\n";
        }
            
        
        outfileHisto.close();*/
    //}

    // eval Ciso and variances
    double Ciso[NTYPE*NS], delta[NTYPE*NS], Delta[NTYPE*NS], Chi4[NTYPE], Chi4Norm[NTYPE];
    for (int type=0; type<NTYPE*NS; type++) {
        Ciso[type] = 0.0;
        delta[type] = 0.0;
        Delta[type] = 0.0;
    }
    for (int type=0; type<NTYPE; type++) {
        Chi4[type] = 0.0;
        Chi4Norm[type] = 0.0;
    }
    
    for (int s = 0; s < NS; s++) {
        for (int i = 0; i < N; i++) {

            // finalized isoconfigurational averages
            dyn_avg[i+s*N] /= (double) NI;
            dyn_avg2[i+s*N] /= (double) NI;

            //double Cval;
            //double Ciso_loc = 0.0;
            //double C2iso_loc = 0.0;
            // Calculate averages (path 2)
            //for (int j = 0; j < NHisto; j++) {
            //    Cval = (j + 0.5)/((double) NHisto);
            //    Ciso_loc += dyn_hist_data[i][j] * Cval;
            //    C2iso_loc += dyn_hist_data[i][j] * Cval * Cval;
                //std::cout << j << " " << dyn_hist_data[i][j] << " ";
            //}
            //std::cout << std::endl;

            // check values
            //std::cout << "AVG direct " << dyn_avg[i] << " AVG indirect " << Ciso_loc << " AVG2 direct " << dyn_avg2[i] << " AVG2 indirect " << C2iso_loc << std::endl;

            // Write into isoconfig storage
            dyn_avg_save[i+s*N][flag*(NT+1)+t] = dyn_avg[i+s*N];
            //if(i+s*N<5) cout <<dyn_avg[i+s*N] << std::endl;

            // calculate delta and Delta
            //std::cout << type_data[i] << std::endl;
            Ciso[type_data[i+s*N]*NS+s] += dyn_avg[i+s*N];
            delta[type_data[i+s*N]*NS+s] += dyn_avg[i+s*N]*dyn_avg[i+s*N];
            Delta[type_data[i+s*N]*NS+s] += dyn_avg2[i+s*N];
        }
    }

    // calc coarse-grained structural descriptors
    int NCG_Dyn=5;
    double mean_den_inherent[NCG_Dyn];
    double mean_rest[NCG_Dyn];
    double Ciso_cg[NTYPE*NCG_Dyn*NS], delta_cg[NTYPE*NCG_Dyn*NS], Delta_cg[NTYPE*NCG_Dyn*NS];
    for (int type=0; type<NTYPE*NCG_Dyn*NS; type++) {
        Ciso_cg[type] = 0.0;
        delta_cg[type] = 0.0;
        Delta_cg[type] = 0.0;
    }
    // cg to be done on the level of the trajectories themselves
    /*for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) { // loop over particles

            for (int c=0; c<NCG_Dyn; c++) {
                mean_den_inherent[c] = 0.0;
                mean_rest[c] = 0.0;   
            }
            //double mean=0.0;
            for (int j=0; j<N;j++) { // loop over particle pairs
                double dr = 0.0, dx;
                for (int d=0; d<NDim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                }
                dr = sqrt(dr);

                for (int c=1; c<NCG_Dyn; c++) {
                    double L = c;
                    double w_inherent = exp(-dr/L);
                    mean_den_inherent[c] += w_inherent;
                    mean_rest[c] += w_inherent*dyn_avg_save[j+s*N][flag*(NT+1)+t];
                }
            }

            for (int c=1; c<NCG; c++) {
                dyn_avg_save_cg[i+s*N][(flag*(NT+1)+t)*NCG_Dyn+c] = mean_rest[c]/mean_den_inherent[c];
            }

            
        }
    }*/

    for (int type=0; type<NTYPE; type++) {
        for (int s = 0; s < NS; s++) {
            double dyn_sys_avg_loc = 0.0;
            //for (int i=0; i<N;i++) if (type==type_data[i+s*N]) dyn_sys_avg_loc += dyn_avg[i+s*N];
            for (int i=0; i<N;i++) 
                if (type==type_data[i+s*N]) {
                    if(DynNames[flag]=="MD") {if(dyn_avg[i+s*N] > 0.438) dyn_sys_avg_loc += 1.0;}
                    else if(DynNames[flag]=="LOG(MD)") {if(dyn_avg[i+s*N] > -0.3) dyn_sys_avg_loc += 1.0;}
                    else if(DynNames[flag]=="LOG(FRES)") {if(dyn_avg[i+s*N] > 1.0) dyn_sys_avg_loc += 1.0;}
                    else if(DynNames[flag]=="LOG(UTH)") {if(dyn_avg[i+s*N] > -0.3) dyn_sys_avg_loc += 1.0;}
                    else {dyn_sys_avg_loc += dyn_avg[i+s*N];}
                }
            Chi4Norm[type]+=dyn_sys_avg_loc;
            Chi4[type]+=dyn_sys_avg_loc*dyn_sys_avg_loc;
        }
        Chi4Norm[type] /= (double) NS;
        Chi4[type] /= (double) NS;
    }

    for (int type=0; type<NTYPE; type++) {
        double R_loc = 0.0;
        double R2_loc = 0.0;
        double Ciso_loc = 0.0;
        double Ciso2_loc = 0.0;
        double Chi4_loc = 0.0;
        double Chi42_loc = 0.0;
        for (int s = 0; s < NS; s++) {
            Ciso[type*NS+s] /= (double) NPerType[type];
            delta[type*NS+s] /= (double) NPerType[type];
            Delta[type*NS+s] /= (double) NPerType[type];
        
            delta[type*NS+s] = delta[type*NS+s] - Ciso[type*NS+s]*Ciso[type*NS+s];
            Delta[type*NS+s] = Delta[type*NS+s] - Ciso[type*NS+s]*Ciso[type*NS+s];

            Ciso_loc += Ciso[type*NS+s];
            Ciso2_loc += Ciso[type*NS+s]*Ciso[type*NS+s];
            R_loc += delta[type*NS+s]/Delta[type*NS+s];
            R2_loc += delta[type*NS+s]/Delta[type*NS+s]*delta[type*NS+s]/Delta[type*NS+s];
       }
        
        dyn_pred[t+NT*flag][5*type + 0] = Ciso_loc/NS;
        dyn_pred[t+NT*flag][5*type + 1] = sqrt(Ciso2_loc/NS - Ciso_loc/NS*Ciso_loc/NS)/sqrt(NS);
        dyn_pred[t+NT*flag][5*type + 2] = R_loc/NS;
        dyn_pred[t+NT*flag][5*type + 3] = sqrt(R2_loc/NS - R_loc/NS*R_loc/NS)/sqrt(NS);
        dyn_pred[t+NT*flag][5*type + 4] = (Chi4[type] - Chi4Norm[type]*Chi4Norm[type])/((double) NPerType[type]);
        //std::cout << Ciso[type] << std::endl;

    }

    // calculate dynamical histograms
    calc_histograms_dynamics(t, flag);

    // evaluate information theory (dynamics)
    //eval_information_theory_dynamics(t);

    // evalue connection dynamics <-> structure
    double * struct_loc;
    for (int j = 0; j < NStructTotal; j++) {
        int cmax = NCG;
        if (j>=struct_read_flag && j < struct_read_flag + struct_read_Ntotal) cmax=1;
        for (int c = 0; c < cmax; c++) {
            struct_loc = struct_local[j*NCG+c];
            // rest histogramms
            for (int j = 0; j < NHisto; j++) {
                for (int type=0; type<NTYPE; type++) {
                    for (int s=0; s<NHistoStruct; s++) {
                        dyn_struct_hist_iso[type][s*NHisto+j] = 0.0;
                        dyn_struct_hist_val[type][s*NHisto+j] = 0.0;
                    }
                }
            }

            // calculate histograms for mutual information (dynamics <-> structure)
            //calc_histograms_information_theory(struct_loc, loc, j, c);

            // evaluate mututal information
            //eval_information_theory_correlation(t,j,c);

            // evaluate pearon and spearman rank coefficients
            eval_pearson_spearman_correlation(t,struct_loc,j,flag,c);

      }
    }

}

void calc_histograms_dynamics(int t, int flag){
    // calc and save histograms for actual motion and isoconfigurational averages
    double hist_lower = dyn_ranges[flag][0];
    double hist_upper = dyn_ranges[flag][1];
    for (int i = 0; i < NS*N; i++) {
        int valint_dyn;
        if (dyn_avg[i] > hist_upper - EPS && dyn_avg[i] < hist_upper + EPS) valint_dyn = NHisto - 1;
        else valint_dyn = (dyn_avg[i]-hist_lower)/(hist_upper - hist_lower)* ((double)NHisto);

        if(valint_dyn >= 0 && valint_dyn < NHisto) dyn_hist_iso[t+(NT+1)*flag][valint_dyn+NHisto*type_data[i]]+= 1.0;
        for (int l = 0; l < NHisto; l++) {
            dyn_hist_val[t+(NT+1)*flag][l+NHisto*type_data[i]] += dyn_hist_data[i][l];
        }
    }

    // normalize histograms
    for (int l = 0; l < NHisto; l++) {
        for (int type=0; type<NTYPE; type++) {
            dyn_hist_val[t+(NT+1)*flag][l+NHisto*type] /= (double) NS*NPerType[type];
            dyn_hist_iso[t+(NT+1)*flag][l+NHisto*type] /= (double) NS*NPerType[type];
        }
    }
}

void eval_information_theory_dynamics(int t) {
    // calc Imot et Is
    /*double I_dyn_struct_rel[NTYPE];
    double I_dyn_struct[NTYPE];
    double I_dyn_tot[NTYPE];
    // then eval the information
    for (int type=0; type<NTYPE; type++) {
        I_dyn_struct_rel[type] = 0.0;
        I_dyn_struct[type] = 0.0;
        I_dyn_tot[type] = 0.0;
    }
    for (int i = 0; i < NS*N; i++) {
        for (int l = 0; l < NHisto; l++) {
            if (dyn_hist_data[i][l] > 0 ) {
                I_dyn_struct_rel[type_data[i]] += dyn_hist_data[i][l]* log2 (dyn_hist_data[i][l]/dyn_hist_val[t][l+NHisto*type_data[i]]);
                I_dyn_struct[type_data[i]] += dyn_hist_data[i][l]* log2 (dyn_hist_data[i][l]);
                
            } 
        }
    }
    for (int type=0; type<NTYPE; type++) {
        for (int l = 0; l < NHisto; l++) {
            if (dyn_hist_iso[t][l+NHisto*type] > 0 ) I_dyn_tot[type] += dyn_hist_iso[t][l+NHisto*type]* log2 (dyn_hist_iso[t][l+NHisto*type]);
        }
    }
    for (int type=0; type<NTYPE; type++) {
        I_dyn_struct_rel[type] /= (double) NS*NPerType[type];
        I_dyn_struct[type] /= (double) NS*NPerType[type];
    }

    
    for (int type=0; type<NTYPE; type++) {
        dyn_pred[t][6*type + 3] = I_dyn_struct_rel[type];
        dyn_pred[t][6*type + 4] = I_dyn_struct[type];
        dyn_pred[t][6*type + 5] = I_dyn_tot[type];
    }*/
}

void calc_histograms_information_theory(double * struct_array, int loc_dyn,int loc_struct ,int c){
    /*double hist_lower_dyn = dyn_ranges[loc_dyn][0];
    double hist_upper_dyn = dyn_ranges[loc_dyn][1];
    double hist_lower_struct = struct_ranges[loc_struct*NCG+c][0];
    double hist_upper_struct = struct_ranges[loc_struct*NCG+c][1];
    for (int i = 0; i < NS*N; i++) {
        int valint_dyn;
        if (dyn_avg[i] > hist_upper_dyn - EPS) valint_dyn = NHisto - 1;
        else valint_dyn = (dyn_avg[i]-hist_lower_dyn)/(hist_upper_dyn - hist_lower_dyn)* ((double)NHisto);
        int valint_struct;
        if (struct_array[i] > hist_upper_struct - EPS) valint_struct = NHistoStruct - 1;
        else valint_struct = (struct_array[i]-hist_lower_struct)/(hist_upper_struct - hist_lower_struct)* ((double)NHistoStruct);
        dyn_struct_hist_iso[type_data[i]][valint_struct*NHisto+valint_dyn]+=1.0;
        for (int l = 0; l < NHisto; l++) {
            dyn_struct_hist_val[type_data[i]][l+NHisto*valint_struct] += dyn_hist_data[i][l];
        }
    }

    // normalize histograms
    for (int type=0; type<NTYPE; type++) {
        for (int l = 0; l < NHisto; l++) {
            for (int s = 0; s < NHistoStruct; s++) {
                dyn_struct_hist_val[type][l+NHisto*s] /= (double) NS*NPerType[type];
                dyn_struct_hist_iso[type][l+NHisto*s] /= (double) NS*NPerType[type];
            }
        }
    }*/
}

void eval_information_theory_correlation(int t,int loc_struct,int c){
    /*for (int type=0; type<NTYPE; type++) {
        for (int l = 0; l < NHisto; l++) {
            for (int s = 0; s < NHistoStruct; s++) {
                //if (type==0 && loc_struct==0 && c==2) std::cout << l << " " << s << " " << dyn_struct_hist_iso[type][l+NHisto*s] << " " << dyn_hist_iso[t][l+NHisto*type] << " " << struct_hist[loc_struct*NCG+c][type*NHistoStruct+s] << std::endl;
                // first Imut for the particle motion
                if (dyn_struct_hist_val[type][l+NHisto*s] > 0 && dyn_hist_val[t][l+NHisto*type] > 0 && struct_hist[loc_struct*NCG+c][type*NHistoStruct+s] > 0)
                    dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal][type]+= dyn_struct_hist_val[type][l+NHisto*s]*log2(dyn_struct_hist_val[type][l+NHisto*s]/(dyn_hist_val[t][l+NHisto*type]*struct_hist[loc_struct*NCG+c][type*NHistoStruct+s]));
                // then Imut for the propernsities
                if (dyn_struct_hist_iso[type][l+NHisto*s] > 0 && dyn_hist_iso[t][l+NHisto*type] > 0 && struct_hist[loc_struct*NCG+c][type*NHistoStruct+s] > 0)
                    dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal][type+NTYPE]+= dyn_struct_hist_iso[type][l+NHisto*s]*log2(dyn_struct_hist_iso[type][l+NHisto*s]/(dyn_hist_iso[t][l+NHisto*type]*struct_hist[loc_struct*NCG+c][type*NHistoStruct+s]));
            }
        }
        dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal][type+4*NTYPE] = dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal][type+NTYPE]/dyn_pred[t][6*type + 3];
        dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal][type+5*NTYPE] = dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal][type+NTYPE]/(dyn_pred[t][6*type + 3]-dyn_pred[t][6*type + 4]);
    }*/
}

void eval_pearson_spearman_correlation(int t,double * struct_array,int loc_struct, int flag, int c){
    // first pearson correlation coefficient
    double mean[2*NTYPE];
    double cov[NTYPE];
    double var_dyn[NTYPE];
    double var_struct[NTYPE];
    for (int type=0; type<NTYPE; type++) {
        mean[type]=0.0;
        mean[type+NTYPE]=0.0;
        cov[type]=0.0;
        var_dyn[type]=0.0;
        var_struct[type]=0.0;
    }
    for (int i = 0; i < NS*N; i++) {
        mean[type_data[i]] += dyn_avg[i];
        mean[type_data[i]+NTYPE] += struct_array[i];
    }
    for (int type=0; type<NTYPE; type++) {
        mean[type] /= (double) NS*NPerType[type];
        mean[type+NTYPE] /= (double) NS*NPerType[type];
    }
    // now calculate covariances and variances
    for (int i = 0; i < NS*N; i++) {
        cov[type_data[i]] += (dyn_avg[i]-mean[type_data[i]])*(struct_array[i]-mean[type_data[i]+NTYPE]);
        var_dyn[type_data[i]] += (dyn_avg[i]-mean[type_data[i]])*(dyn_avg[i]-mean[type_data[i]]);
        var_struct[type_data[i]] += (struct_array[i]-mean[type_data[i]+NTYPE])*(struct_array[i]-mean[type_data[i]+NTYPE]);
    }
    for (int type=0; type<NTYPE; type++) {
        dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal+NCG*NStructTotal*NT*flag][type]= cov[type]/sqrt(var_dyn[type]*var_struct[type]);
    }

    // calculate ranks
    double save_dyn[NS*N]; 
    double save_struct[NS*N]; 
    for(int i = 0; i < NS*N; i++) {
        save_dyn[i] = dyn_avg[i]; 
        save_struct[i] = struct_array[i]; 
    }
    sort(save_dyn, save_dyn+NS*N); // sorting the array 
    sort(save_struct, save_struct+NS*N); // sorting the array 
    map<double, int> rank_dyn; 
    map<double, int> rank_struct; 
    for (int i = 0; i < NS*N; i++) { 
        rank_dyn[save_dyn[i]] = i; 
        rank_struct[save_struct[i]] = i; 
    } 
    /*for (int i = 0; i < NS*N; i++) {
        if (rank_struct[struct_array[i]]==0) cout << struct_array[i] << " "; 
        if (rank_struct[struct_array[i]]==1) cout << struct_array[i] << " "; 
    }
    cout << endl;*/

    // now spearman
    for (int type=0; type<NTYPE; type++) {
        mean[type]=0.0;
        mean[type+NTYPE]=0.0;
        cov[type]=0.0;
        var_dyn[type]=0.0;
        var_struct[type]=0.0;
    }
    for (int i = 0; i < NS*N; i++) {
        mean[type_data[i]] += rank_dyn[dyn_avg[i]];
        mean[type_data[i]+NTYPE] += rank_struct[struct_array[i]];
    }
    for (int type=0; type<NTYPE; type++) {
        mean[type] /= (double) NS*NPerType[type];
        mean[type+NTYPE] /= (double) NS*NPerType[type];
    }
    // now calculate covariances and variances
    for (int i = 0; i < NS*N; i++) {
        cov[type_data[i]] += (rank_dyn[dyn_avg[i]]-mean[type_data[i]])*(rank_struct[struct_array[i]]-mean[type_data[i]+NTYPE]);
        var_dyn[type_data[i]] += (rank_dyn[dyn_avg[i]]-mean[type_data[i]])*(rank_dyn[dyn_avg[i]]-mean[type_data[i]]);
        var_struct[type_data[i]] += (rank_struct[struct_array[i]]-mean[type_data[i]+NTYPE])*(rank_struct[struct_array[i]]-mean[type_data[i]+NTYPE]);
    }
    for (int type=0; type<NTYPE; type++) {
        dyn_struct_pred[loc_struct*NCG+c+t*NCG*NStructTotal+NCG*NStructTotal*NT*flag][type+NTYPE]= cov[type]/sqrt(var_dyn[type]*var_struct[type]);
    }
}


// calculate reaarangement timescale from isoconf average
void eval_timescale(int flag, double threshold){

    for (int i = 0; i < NS*N; i++) {
        dyn_avg_save[i][flag*(NT+1)+NT] = 0.0;

        for (int t=1; t<NT; t++) {
            //cout << t << " " <<  dyn_avg_save[i][flag*(NT+1)+t] << std::endl;
            if(dyn_avg_save[i][flag*(NT+1)+t] < threshold && dyn_avg_save[i][flag*(NT+1)+NT] == 0.0) {
                dyn_avg_save[i][flag*(NT+1)+NT] = time_data[t-1]*timestep - (dyn_avg_save[i][flag*(NT+1)+t-1] - threshold)*(time_data[t-1]*timestep-time_data[t]*timestep)/(dyn_avg_save[i][flag*(NT+1)+t-1]-dyn_avg_save[i][flag*(NT+1)+t]);
                //cout << t << " " <<  time_data[t-1]*timestep << " "  << time_data[t]*timestep << " " << dyn_avg_save[i][flag*(NT+1)+NT] << std::endl;
                break;
            }
        }
    }

    // several particles will have a time scale larger than the largest simulated time scale
    // solution: average over all particles with the longest timescales which were still recorded ( NT - 2 < t < NT -1 )
    // -> this gives a master curve which can be used to extrapolate by time-temperature superposition the very slow particles
    double dyn_avg_super[NT];
    for (int t=1; t<NT; t++) {
        dyn_avg_super[t] = 0.0;
    }
    int count=0;
    int tmax = NT-1;
    //cout << time_data[NT-2] << " " << time_data[NT-1] << " " << tmax << " " << NT << std::endl;
    for (int i = 0; i < NS*N; i++) {
        if (dyn_avg_save[i][flag*(NT+1)+NT] > 0.95*time_data[tmax]*timestep ) {
            count++;
            for (int t=1; t<NT; t++) {
                dyn_avg_super[t] += dyn_avg_save[i][flag*(NT+1)+t];
            }
        }
    }
    for (int t=1; t<NT; t++) {
        dyn_avg_super[t] /= (double) count;
        //cout << t << " " << time_data[t]*timestep << " "<< dyn_avg_super[t] << std::endl;
    }

    for (int i = 0; i < NS*N; i++) {
        // find time scale for all remaining particles
        if (dyn_avg_save[i][flag*(NT+1)+NT] < 0.0001) {
            double val = dyn_avg_save[i][flag*(NT+1)+NT-1];
            if (val > 0.98) val = 0.98;  //introduce cutoff to avoid articifically long tails
            for (int t=1; t<NT; t++) {
                if(val > dyn_avg_super[t]) {
                    dyn_avg_save[i][flag*(NT+1)+NT] = time_data[t-1]*timestep - (dyn_avg_super[t-1] - val)*(time_data[t-1]*timestep-time_data[t]*timestep)/(dyn_avg_super[t-1]-dyn_avg_super[t]);
                    //cout << t << " " <<  time_data[t-1]*timestep << " "  << time_data[t]*timestep << " " << dyn_avg_save[i][flag*(NT+1)+NT] << std::endl;
                    break;
                }
            }
            dyn_avg_save[i][flag*(NT+1)+NT] = time_data[NT-1]*timestep*time_data[NT-1]*timestep/dyn_avg_save[i][flag*(NT+1)+NT];
            
            //cout << dyn_avg_save[i][flag*(NT+1)+NT] << " " << val << std::endl; 
        }
    }

}

void print_isoconf(int flag){

    QString pathOrig = QString::fromStdString(folderOut);

    // print predictabilities
    QString pathPred = pathOrig;
    pathPred.append(QString("/isoconf_predictability_%1.dat").arg(QString::fromStdString(DynNames[flag])));
    QFile outfilePred(pathPred);   // input file with xyz
    outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred(&outfilePred);
    for (int t=1; t<NT; t++) {
        outPred << time_data[t]*timestep << " ";
        for (int type=0; type<NTYPE; type++) {
            for (int j = 0; j < 5; j++) {
                outPred << dyn_pred[t+NT*flag][5*type + j] << " ";
            }
        }
        outPred << "\n";
    }
    outfilePred.close();

    // print information theory and correlation coefficients
    for (int j = 0; j < NStructTotal; j++) {
        int cmax = NCG;
        if (j>=struct_read_flag && j < struct_read_flag + struct_read_Ntotal) cmax=1;
        for (int c = 0; c < cmax; c++) {
            QString pathPred = pathOrig;
            pathPred.append(QString("/isoconf_correlation_%1_%2_%3.dat").arg(QString::fromStdString(DynNames[flag])).arg(QString::fromStdString(StructNames[j])).arg(c));
            QFile outfilePred(pathPred);   // input file with xyz
            outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
            QTextStream outPred(&outfilePred);
            outPred << "Time Kp Ks Rs Rtot\n";
            for (int t=1; t<NT; t++) {
                outPred << time_data[t]*timestep << " ";
                for (int type=0; type<NTYPE; type++) {
                    for (int k = 0; k < 1; k++) { // only print pearson
                        outPred << dyn_struct_pred[j*NCG+c+t*NCG*NStructTotal+NCG*NStructTotal*NT*flag][type+k*NTYPE] << " " << dyn_struct_pred[j*NCG+c+t*NCG*NStructTotal+NCG*NStructTotal*NT*flag][type+k*NTYPE]*dyn_pred[t+NT*flag][4*type + 2] << " ";
                    }
                }
                outPred << "\n";
            }
            outfilePred.close();
        }
    }


    // print histograms
    double hist_lower = dyn_ranges[flag][0];
    double hist_upper = dyn_ranges[flag][1];
    QString pathHisto = pathOrig;
    pathHisto.append(QString("/isoconf_histograms_%1.dat").arg(QString::fromStdString(DynNames[flag])));
    QFile outfileHisto(pathHisto);   // input file with xyz
    outfileHisto.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outHisto(&outfileHisto);
    for (int t=1; t<=NT; t++) {
        if (t<NT) outHisto << time_data[t]*timestep << "\n";
        else outHisto << "Relaxation time" << "\n";
        for (int j = 0; j < NHisto; j++) {
            if (t<NT) outHisto << hist_lower + (j+0.5)/((double)NHisto)*(hist_upper - hist_lower) << " ";
            else outHisto << hist_lower_time + (j+0.5)/((double)NHisto)*(hist_upper_time - hist_lower_time) << " ";
            for (int type=0; type<NTYPE; type++) {
                outHisto <<  dyn_hist_val[t+(NT+1)*flag][j+NHisto*type] << " ";
            }
            for (int type=0; type<NTYPE; type++) {
                outHisto <<  dyn_hist_iso[t+(NT+1)*flag][j+NHisto*type] << " ";
            }
            outHisto << "\n";
        }
        outHisto << "\n\n";
    }
    outfileHisto.close();


}

void print_traj(double * save_dyn, int flag){

    QString pathOrig = QString::fromStdString(folderOut);

    // print only first 8 files
    int NIloc = NI;
    if (NIloc > 8) NIloc = 8;
    for (int j=0; j<NIloc; j++) {

        QString pathPred = pathOrig;
        pathPred.append(QString("struct_trajectory_%1_%2.xyz").arg(j).arg(QString::fromStdString(DynNames[flag])));
        std::cout << pathPred.toStdString() << std::endl;
        QFile outfilePred(pathPred);   // input file with xyz
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred(&outfilePred);
        for (int t=0; t<NT; t++) {

            // calc mean
            /*double mean = 0.0;
            for (int i = 0; i < N; i++) {
                mean+=save_dyn[t*N*NI+j*N+i];
            }
            mean /= (double) N;
            for (int i = 0; i < N; i++) {
                if (save_dyn[t*N*NI+j*N+i] < mean) save_dyn[t*N*NI+j*N+i] = 0;
                else save_dyn[t*N*NI+j*N+i] = 1;
            }*/

            outPred << N << "\n";
            outPred << "Properties=species:I:1:pos:R:" << NDim;
            outPred << ":LOG(FRES)/BB:R:1";
            outPred << " time " << time_data[t]*timestep << "\n";
            for (int i = 0; i < N; i++) {
                if (NDim == 2) {
                    outPred << type_data[i]+1 << " " << xyz_data[i][t*NDim+NDim*NT*j] << " " << xyz_data[i][1+t*NDim+NDim*NT*j] << " ";
                } else {
                    outPred << type_data[i]+1 << " " << xyz_data[i][t*NDim+NDim*NT*j] << " " << xyz_data[i][1+t*NDim+NDim*NT*j] << " " << xyz_data[i][2+t*NDim+NDim*NT*j] << " ";
                }
                outPred << save_dyn[t*N*NI+j*N+i] << " ";
                outPred << "\n";
            }
        }
        outfilePred.close();
    }

}


void print_R(double * save_dyn,int flag){

    double * Cg_avg = new double [NT*NTYPE];
    double * Cg_iso2 = new double [NT*NTYPE];
    double * Cg2_iso = new double [NT*NTYPE];
    double * chi4_iso = new double [NT*NTYPE];
    double * chi4 = new double [NT*NTYPE];
    double * R4 = new double [NT*NTYPE];

    QString pathOrig = QString::fromStdString(folderOut);
    QString pathPred2 = pathOrig;
    pathPred2.append(QString("isoconf_trajectories.dat"));
    QFile outfilePred2(pathPred2);   // input file with xyz
    outfilePred2.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred2(&outfilePred2);

    for (int type=0; type<NTYPE; type++) {
        for (int t=1; t<NT; t++) {
            if (type==0) outPred2 << time_data[t]*timestep << " ";
            Cg_avg[t+NT*type] = 0.0;
            Cg_iso2[t+NT*type] = 0.0;
            Cg2_iso[t+NT*type] = 0.0;
            chi4_iso[t+NT*type] = 0.0;
            chi4[t+NT*type] = 0.0;
            R4[t+NT*type] = 0.0;
            for (int s=0; s<NS; s++) {
                double Cg_iso2_loc = 0.0;
                for (int j=0; j<NI; j++) {
                    double Cg = 0.0;
                    for (int i=0; i<N; i++) {
                        if(type_data[i+s*N] == type) Cg += save_dyn[s*NT*NI*N+t*N*NI+j*N+i];
                    }
                    Cg /= (double) NPerType[type];
                    if (type==0) outPred2 << Cg << " ";
                    Cg_iso2_loc += Cg;
                    Cg_avg[t+NT*type] += Cg;
                    Cg2_iso[t+NT*type] += Cg*Cg;
                }
                Cg_iso2[t+NT*type] += Cg_iso2_loc*Cg_iso2_loc / ( (double) NI*NI );
            }
            Cg_iso2[t+NT*type] /= (double) NS;
            Cg_avg[t+NT*type] /= (double) NS*NI;
            Cg2_iso[t+NT*type] /= (double) NS*NI;
            chi4[t+NT*type] = NPerType[type]*(Cg2_iso[t+NT*type] - Cg_avg[t+NT*type]*Cg_avg[t+NT*type]);
            chi4_iso[t+NT*type] = NPerType[type]*(Cg_iso2[t+NT*type] - Cg_avg[t+NT*type]*Cg_avg[t+NT*type]);
            if (t==36) printf("BB/chi4/var/mean %f %f %f %f\n",Cg_avg[t+NT*type],chi4_iso[t+NT*type],Cg_iso2[t+NT*type]*NPerType[type]*NPerType[type], Cg_avg[t+NT*type]*NPerType[type] );
            R4[t+NT*type] = chi4_iso[t+NT*type] / chi4[t+NT*type];
            if (type==0) outPred2 << "\n";
        }
    }
    outfilePred2.close();

    // print R factor
    QString pathPred = pathOrig;
    pathPred.append(QString("/isoconf_R_%1.dat").arg(QString::fromStdString(DynNames[flag])));
    QFile outfilePred(pathPred);   // input file with xyz
    outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred(&outfilePred);
    for (int t=1; t<NT; t++) {
        outPred << time_data[t]*timestep << " ";
        for (int type=0; type<NTYPE; type++) {
            outPred << dyn_pred[t+NT*flag][5*type + 0] << " ";
            outPred << dyn_pred[t+NT*flag][5*type + 2] << " ";
            outPred << chi4[t+NT*type] << " ";
            outPred << chi4_iso[t+NT*type] << " ";
            outPred << R4[t+NT*type] << " ";
        }
        outPred << "\n";
    }
    outfilePred.close();
}