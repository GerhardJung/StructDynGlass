 
#include "eval_isoconf.h"
#include "defs.h"

void reset_dyn_hist(){
    for (int i = 0; i < NS*N; i++) {
        for (int j = 0; j < NHisto; j++) {
            dyn_hist_data[i][j] = 0.0;
        }
    }
}

void add_histogram(int s, int i, double val){
    int valint = val*NHisto;
    dyn_hist_data[s*N+i][valint]++;
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

void eval_isoconf(int t){

    norm_histogram();

    //calc Imot et Is
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
                Imot[type_data[i]] -= dyn_hist_data[i][j]* log (den[j]);
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
        dyn_pred[t][5*type + 3] = Is[type];
    }
}