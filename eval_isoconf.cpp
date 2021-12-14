 
#include "eval_isoconf.h"
#include "defs.h"

#define EPS 0.000000001

void reset_dyn(){
    for (int i = 0; i < NS*N; i++) {
        for (int j = 0; j < NHisto; j++) {
            dyn_hist_data[i][j] = 0.0;
        }
        dyn_avg[i] = 0.0;
        dyn_avg2[i] = 0.0;
    }
}

void add_histogram_avg(int s, int i, double hist_lower, double hist_upper, double val){
    int valint;
    if (val > hist_upper - EPS) valint = NHisto - 1;
    else valint = (val-hist_lower)/(hist_upper - hist_lower)* ((double)NHisto);
    //if ( valint < 0 || valint >= NHisto) std::cout << valint << " " << val << " " <<  (val-hist_lower)/(hist_upper - hist_lower) << " " <<  hist_upper << " " <<  hist_lower << std::endl;
    //std::cout << val << " " << valint << std::endl;
    dyn_hist_data[s*N+i][valint]++;

    dyn_avg[s*N+i] += val;
    dyn_avg2[s*N+i] += val*val; 
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

void eval_isoconf(int t, std::string dyn){

    norm_histogram();

    // eval Ciso and variances
    double Ciso[NTYPE], delta[NTYPE], Delta[NTYPE];
    for (int type=0; type<NTYPE; type++) {
        Ciso[type] = 0.0;
        delta[type] = 0.0;
        Delta[type] = 0.0;
    }

    for (int i = 0; i < NS*N; i++) {

        // finalized isoconfigurational averages
        dyn_avg[i] /= (double) NI;
        dyn_avg2[i] /= (double) NI;

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
        if (dyn=="BB") dyn_bb_avg[i][t] = dyn_avg[i];
        if (dyn=="EXP") dyn_exp_avg[i][t] = dyn_avg[i];
        if (dyn=="ISFX") dyn_isf_avg[i][t*dim] = dyn_avg[i];
        if (dyn=="ISFY") dyn_isf_avg[i][t*dim+1] = dyn_avg[i];
        if (dyn=="ISFZ") dyn_isf_avg[i][t*dim+2] = dyn_avg[i];

        // calculate delta and Delta
        //std::cout << type_data[i] << std::endl;
        Ciso[type_data[i]] += dyn_avg[i];
        delta[type_data[i]] += dyn_avg[i]*dyn_avg[i];
        Delta[type_data[i]] += dyn_avg2[i];
    }

    for (int type=0; type<NTYPE; type++) {
        Ciso[type] /= (double) NS*NPerType[type];
        delta[type] /= (double) NS*NPerType[type];
        Delta[type] /= (double) NS*NPerType[type];
        //if (type ==0) std::cout << Ciso[type] << " " << delta[type] << " " << Delta[type] << " " << std::endl; 
        delta[type] = delta[type] - Ciso[type]*Ciso[type];
        Delta[type] = Delta[type] - Ciso[type]*Ciso[type];
        dyn_pred[t][5*type + 0] = Ciso[type];
        dyn_pred[t][5*type + 1] = delta[type];
        dyn_pred[t][5*type + 2] = Delta[type];
        //std::cout << Ciso[type] << std::endl;
    }


    // calc Imot et Is
    double den[NHisto*NTYPE];
    double Imot[NTYPE];
    double Is[NTYPE];
    // first calculate denominator
    for (int j = 0; j < NHisto; j++) {
        for (int type=0; type<NTYPE; type++) {
            den[j+NHisto*type] = 0.0;
        }
        
        for (int i = 0; i < NS*N; i++) {
            den[j+NHisto*type_data[i]] += dyn_hist_data[i][j];
        }
        for (int type=0; type<NTYPE; type++) {
            den[j+NHisto*type] /= (double) NS*NPerType[type];
        }
    }
    // then eval the information
    for (int type=0; type<NTYPE; type++) {
        Imot[type] = 0.0;
        Is[type] = 0.0;
    }
    for (int i = 0; i < NS*N; i++) {
        for (int j = 0; j < NHisto; j++) {
            if (dyn_hist_data[i][j] > 0 ) {
                Imot[type_data[i]] -= dyn_hist_data[i][j]* log (den[j+NHisto*type_data[i]]);
                Is[type_data[i]] += dyn_hist_data[i][j]* log (dyn_hist_data[i][j]);
            } 
        }
    }
    for (int type=0; type<NTYPE; type++) {
        Imot[type] /= (double) NS*NPerType[type];
        Is[type] /= (double) NS*NPerType[type];
        Is[type] = Imot[type] + Is[type];
    }

    
    for (int type=0; type<NTYPE; type++) {
        dyn_pred[t][5*type + 3] = Imot[type];
        dyn_pred[t][5*type + 4] = Is[type];
    }

    // calc and save mean histograms
    for (int i = 0; i < NS*N; i++) {
        int valint = dyn_avg[i]*NHisto;
        dyn_hist_iso[t][valint+NHisto*type_data[i]]++;
    }
    for (int j = 0; j < NHisto; j++) {
        for (int type=0; type<NTYPE; type++) {
            dyn_hist_val[t][j+NHisto*type] = den[j+NHisto*type];
            dyn_hist_iso[t][j+NHisto*type] /= (double) NS*NPerType[type];
        }
    }
}

void print_isoconf(std::string dyn){

    QString pathOrig = QString::fromStdString(folderOut);

    // print predictabilities
    QString pathPred = pathOrig;
    pathPred.append(QString("/isoconf_predictability_%1.dat").arg(QString::fromStdString(dyn)));
    QFile outfilePred(pathPred);   // input file with xyz
    outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outPred(&outfilePred);
    for (int t=1; t<NT; t++) {
        outPred << time_data[t]*timestep << " ";
        for (int type=0; type<NTYPE; type++) {
            for (int j = 0; j < 5; j++) {
                outPred << dyn_pred[t][5*type + j] << " ";
            }
        }
        outPred << "\n";
    }
    outfilePred.close();

    // print histograms
    QString pathHisto = pathOrig;
    pathHisto.append(QString("/isoconf_histograms_%1.dat").arg(QString::fromStdString(dyn)));
    QFile outfileHisto(pathHisto);   // input file with xyz
    outfileHisto.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outHisto(&outfileHisto);
    for (int t=1; t<NT; t++) {
        outHisto << time_data[t]*timestep << "\n";
        for (int j = 0; j < NHisto; j++) {
            if(dyn == "BB") outHisto << bb_hist_lower + (j+0.5)/((double)NHisto)/(bb_hist_upper - bb_hist_lower) << " ";
            if(dyn == "EXP") outHisto << exp_hist_lower + (j+0.5)*NHisto/(exp_hist_upper - exp_hist_lower) << " ";
            if(dyn == "ISFX" || dyn == "ISFY" || dyn == "ISFZ" ) outHisto << isf_hist_lower + (j+0.5)*NHisto/(isf_hist_upper - isf_hist_lower) << " ";
            for (int type=0; type<NTYPE; type++) {
                outHisto <<  dyn_hist_val[t][j+NHisto*type] << " ";
            }
            for (int type=0; type<NTYPE; type++) {
                outHisto <<  dyn_hist_iso[t][j+NHisto*type] << " ";
            }
            outHisto << "\n";
        }
        outHisto << "\n\n";
    }
    outfileHisto.close();

}