
#include "struct_epot.h"
#include "defs.h"
#include "pbc.h"

void eval_struct_epot(){

    std::cout << "EVAL STRUCT EPOT " << std::endl; 

    int iType,jType;
    double epot = 0.0;
    double epot_inh = 0.0;

    for (int s=0; s<NS;s++) { // loop over structures

        for (int i=0; i<N;i++) { // loop over particles
            int ncount = 0;
            iType = type_data[i+s*N];

            for (int j=0; j<N;j++) { // loop over particle pairs
                jType = type_data[j+s*N];
                double dr = 0.0, dx;
                double dr_inherent = 0.0, dx_inherent;
                for (int d=0; d<NDim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                    dx_inherent = xyz_inherent_data[i+s*N][d] - xyz_inherent_data[j+s*N][d];
                    apply_pbc(dx_inherent);
                    dr_inherent += dx_inherent*dx_inherent;
                }

                // calculate epot
                if (i!=j) epot += 0.5*calc_epot2(i+s*N, j+s*N, dr);
                if (i!=j) epot_inh += 0.5*calc_epot2(i+s*N, j+s*N, dr_inherent);
                
            } 

        }
    }

	epot /= (double) NS * N;
	epot_inh /= (double) NS * N;
	
	    std::cout << "EPOT=" << epot << " EPOT(INH)=" << epot_inh << std::endl; 

}


// help function
double determine_sigma2(int i, int j) {
    if(model=="KA2-2D") {
        int iType = type_data[i];
        int jType = type_data[j];
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.88;
        if (iType==2 && jType == 2) return 0.94;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 0.8;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.9;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 0.8;
    }
    if (model=="KA2" || model=="KA" || model=="KASHIBA") {
        int iType = type_data[i];
        int jType = type_data[j];
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.88;
        if (iType==2 && jType == 2) return 0.94;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 0.8;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.9;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 0.84;
    }
    if (model=="POLY") {
        double iRadius = dia_data[i];
        double jRadius = dia_data[j];
        return (iRadius+jRadius)/2.0*(1.0-0.2*abs(iRadius-jRadius));
    }
    return 0.0;
}
double determine_epsilon2(int iType, int jType) {
    if(model=="KA2-2D") {
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.5;
        if (iType==2 && jType == 2) return 0.75;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 1.5;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 0.75;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 1.5;
    } 
    if(model=="KA2" || model=="KA" || model=="KASHIBA") {
        if (iType==0 && jType == 0) return 1.0;
        if (iType==1 && jType == 1) return 0.5;
        if (iType==2 && jType == 2) return 0.75;
        if ((iType==1 && jType == 0) || (iType==0 && jType == 1)) return 1.5;
        if ((iType==2 && jType == 0) || (iType==0 && jType == 2)) return 1.25;
        if ((iType==2 && jType == 1) || (iType==1 && jType == 2)) return 1.0;
    } 
    if(model=="POLY") {
        return 1.0;
    } 
    return 0.0;
}

double calc_epot2(int i, int j, double dist2) {
    double sigma = determine_sigma2(i, j);
    double epsilon = determine_epsilon2(type_data[i], type_data[j]);
    double sigma2=sigma*sigma;

    if(model=="KA2" || model=="KA2-2D") {
        if (dist2 < RC2LJ*sigma2 ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(C0LJ+C2LJ*rij2+C4LJ*rij4 - rij6 + rij12 ) ;
        } 
    } if(model=="KA") {
        if (dist2 < RC2LJ ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(- rij6 + rij12 ) ;
        } 
    } if(model=="KASHIBA") {
        if (dist2 < RC2LJ*sigma2 ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(- rij6 + rij12 -c0SHIBA - ( sqrt(dist2)/sigma - 2.5)*c1SHIBA ) ;
        } 
    } if(model=="POLY") {
        if (dist2 < RC2POLY*sigma2 ) {
            double rij2 = dist2/(sigma2);
            double rij4 = rij2*rij2;
            double rij6 = 1.0/(rij4*rij2);
            double rij12 = rij6*rij6;
            return 4.0*epsilon*(C0POLY+C2POLY*rij2+C4POLY*rij4 + rij12 ) ;
        } 
    } else return 0.0;
}
