
#include "eval_struct.h"
#include "struct_ml.h"
#include "pbc.h"
#include "defs.h"

void eval_struct(double * input_tensor,std::string input_name, int first){
    double Chi4[NTYPE]={}, Chi4Norm[NTYPE]={};
    std::string loc=input_name;
    //loc.resize(loc.size() - 2);

    // calc chi4
    for (int type=0; type<NTYPE; type++) {
        for (int s = 0; s < NS; s++) {
            double struct_sys_avg = 0.0;
            for (int i=0; i<N;i++)
                if (type==type_data[i+s*N]) {
                    if(loc.substr( loc.length() - 2 )=="MD") {if(input_tensor[i+s*N] > 0.438) struct_sys_avg += 1.0;}
                    else if(loc.substr( loc.length() - 2 )=="BB") { struct_sys_avg += 1.0 - input_tensor[i+s*N];}
                    else if(loc.length() > 7 && loc.substr( loc.length() - 7 )=="LOG(MD)") {if(input_tensor[i+s*N] > -0.3) struct_sys_avg += 1.0;}
                    else if(loc.length() > 9 && loc.substr( loc.length() - 9 )=="LOG(FRES)") {if(input_tensor[i+s*N] > 1.0) struct_sys_avg += 1.0;}
                    else if(loc.length() > 8 && loc.substr( loc.length() - 8 )=="LOG(UTH)") {if(input_tensor[i+s*N] > -0.3) struct_sys_avg += 1.0;}
                    else {struct_sys_avg += input_tensor[i+s*N];}
                }

            Chi4Norm[type]+=struct_sys_avg;
            Chi4[type]+=struct_sys_avg*struct_sys_avg;
            //std::cout << StructNames[flag] << " " << struct_sys_avg << std::endl;
        }
        Chi4Norm[type] /= (double) NS;
        Chi4[type] /= (double) NS;
    }

    std::cout << "EVAL STRUCT READ: FINISHED CHI4 " << std::endl; 

    // print chi4
    QString pathPred = QString::fromStdString(folderOut);
    /*std::cout << StructNames[flag].substr( StructNames[flag].length() - 2 ) << std::endl; 
    int t=stoi(StructNames[flag].substr( StructNames[flag].length() - 2 ));

    pathPred.append(QString("/structure_%1.dat").arg(QString::fromStdString(loc)));
    if (t==11) { // delete file
        QFile outfilePred(pathPred);  
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        outfilePred.close();
    }

    std::cout << "EVAL STRUCT READ: FINISHED CHI4 2 " << std::endl; 

    QFile outfilePred2(pathPred);   
    outfilePred2.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append);
    QTextStream outPred(&outfilePred2);
    outPred << time_data[t]*timestep << " ";
    for (int type=0; type<NTYPE; type++) {
            outPred << (Chi4[type] - Chi4Norm[type]*Chi4Norm[type])/((double) NPerType[type]) << " ";
    }
    outPred << "\n";
    
    outfilePred2.close();*/

    double res = 0.0;
    int countres = 0;
    for (int i=0; i<N*NS;i++) {
        if (type_data[i]==0) {
          res += input_tensor[i];
          countres += 1;
        }
    }
    printf("mean: %f\n",res/((double) countres));


    std::cout << "EVAL STRUCT READ: S4 " << std::endl; 

    // calc S4
		double qmin = 2.0 * 3.1415926 / boxL;
		int Nq = 15;
				
		// calc density profiles
		double ** resVec;
    resVec = dmatrix(0,NS*N-1,0,4*Nq-1);
    double c,s,c0,c1,s1,c2,s2;
    for (int i=0; i<N*NS;i++) {
        if (type_data[i]==0) {
          for (int d=0; d<2;d++) {
              for (int m=0; m<Nq;m++) {
                
                if (m == 0) {
                  c = cos(qmin*xyz_data[i][d]);
                  s = sin(qmin*xyz_data[i][d]);
                  //if (i==1 ) {
                  //  printf("%f\n",qmin*xyz_data[i][d]);
                  //}
                  c0 = c;
                }
                else if (m == 1) {
                  c1 = c;
                  s1 = s;
                  c = 2.0 * c0*c1 - 1.0;
                  s = 2.0 * c0*s1;
                } else {
                  c2 = c1;
                  s2 = s1;
                  c1 = c;
                  s1 = s;
                  c = 2.0 * c0 * c1 - c2;
                  s = 2.0 * c0 * s1 - s2;
                }

                resVec[i][2*d*Nq+2*m]= c;
                resVec[i][2*d*Nq+2*m+1]= s;
              }
            } 
        }
    }

    printf("%f\n",resVec[1][10]);
				
    double Re1_pred[(Nq+1)*(Nq+1)*NS]={};    
    double Re2_pred[(Nq+1)*(Nq+1)*NS]={};    
    double Im1_pred[(Nq+1)*(Nq+1)*NS]={};    
    double Im2_pred[(Nq+1)*(Nq+1)*NS]={};    
				
    double Rex,Imx,Rey,Imy;
		for (int Nx=0; Nx<Nq+1;Nx++) {
				for (int Ny=0; Ny<Nq+1;Ny++) {
					 for (int s=0; s<NS;s++) {
              for (int i=0; i<N;i++) {
                if (type_data[i+s*N]==0) {
                  if (Nx == 0) {
                    Rex = 1.0;
                    Imx = 0.0;
                  } else {
                    Rex = resVec[s*N+i][ 2*(Nx-1)];
                    Imx = resVec[s*N+i][ 2*(Nx-1) + 1];
                  }
                  if (Ny == 0) {
                    Rey = 1.0;
                    Imy = 0.0;
                  } else {
                    Rey = resVec[s*N+i][ 2*Nq + 2*(Ny-1)];
                    Imy = resVec[s*N+i][ 2*Nq + 2*(Ny-1) + 1];
                  }

                  if (s==0 && Nx==1 && Ny==2 && (i<15) ) {
                    printf("%f %f %f\n",Rex,Rey,input_tensor[i+s*N]);
                  }

                  Re1_pred[(Nx)*(Nq+1)*NS+(Ny)*NS+s] += (input_tensor[i+s*N])*Rex*Rey;
                  Re2_pred[(Nx)*(Nq+1)*NS+(Ny)*NS+s] -= (input_tensor[i+s*N])*Imx*Imy;
                  Im1_pred[(Nx)*(Nq+1)*NS+(Ny)*NS+s] += (input_tensor[i+s*N])*Rex*Imy;
                  Im2_pred[(Nx)*(Nq+1)*NS+(Ny)*NS+s] += (input_tensor[i+s*N])*Imx*Rey;
                }

              }

              if (s==0 && Nx==1 && Ny==2) {
                    printf("%f\n",Re1_pred[(Nx)*(Nq+1)*NS+(Ny)*NS+s]);
                }


           }
        }
    }
    free_dmatrix(resVec,0,NS*N-1,0,4*Nq-1);
				
		// calc structure factor
    double strucFac_pred[(Nq+1)*(Nq+1)]={};    
		for (int Nx=0; Nx<Nq+1;Nx++) {
				for (int Ny=0; Ny<Nq+1;Ny++) {
					 for (int s=0; s<NS;s++) {
								
							double c_pred = Re1_pred[ (Nx*(Nq+1)+Ny)*NS + s] + Re2_pred[ (Nx*(Nq+1)+Ny)*NS + s] ;
							double s_pred = Im1_pred[ (Nx*(Nq+1)+Ny)*NS + s] + Im2_pred[ (Nx*(Nq+1)+Ny)*NS + s] ;
							strucFac_pred[Nx*(Nq+1)+Ny] += c_pred*c_pred + s_pred*s_pred;
           }

                         if ( Nx==1 && Ny==2) {
                    printf("%f\n",strucFac_pred[Nx*(Nq+1)+Ny]);
                }
        }
    }
						
	 // print
    pathPred = QString::fromStdString(folderOut);
    pathPred.append(QString("/S4_%1.dat").arg(QString::fromStdString(input_name)));

    QFile outfilePred3(pathPred);  
    outfilePred3.open(QIODevice::WriteOnly | QIODevice::Text );
    QTextStream outPred3(&outfilePred3);

		for (int Nx=0; Nx<Nq+1;Nx++) {
				for (int Ny=0; Ny<Nq+1;Ny++) {
						if (Nx != 0 or Ny != 0) {
              double qx = (Nx)*qmin;
              double qy = (Ny)*qmin;
              outPred3 <<  sqrt(qx*qx + qy*qy) << " " << strucFac_pred[Nx*(Nq+1)+Ny]/(NPerType[0]*NS);
              outPred3 << "\n";
            }
        }
    }
    outfilePred3.close();


    // g_4
    std::cout << "EVAL STRUCT READ: G4 " << std::endl; 
		int		NBin = 500;
		double		DeltaBin = 0.04;
		double G4[NBin]={};   	
		double sum_pred=0.0;

    // calculate constants for cell list calculation
    if (struct_ml_flag<0 && (first==1) ) {
      rc = NBin*DeltaBin;
      Ncell = (int) boxL/rc;
      Nmax = 800;
      cell_list_index = imatrix(0,NS-1,0,N-1);
      cell_list = imatrix(0,NS*Ncell*Ncell-1,0,Nmax-1);
      create_cell_lists();
    }
				
		// calc distance matrices, norms and true label results	
		for (int s=0; s<NS;s++) {

				// calc distances
        for (int i=0; i<N;i++) { // loop over particles
            if (type_data[i+s*N]==0) {

              int cell_list_indexi = cell_list_index[s][i];
              int cell_list_y  = cell_list_indexi % Ncell;
              int cell_list_x  = (cell_list_indexi - cell_list_y)/Ncell;
              
              for (int a=cell_list_x - 1; a<= cell_list_x + 1; a++) {
                  for (int b=cell_list_y - 1; b<= cell_list_y + 1; b++) {
                      int aloc = a;
                      int bloc = b;
                      if (aloc < 0) aloc += Ncell;
                      if (aloc >= Ncell) aloc -= Ncell;
                      if (bloc < 0) bloc += Ncell;
                      if (bloc >= Ncell) bloc -= Ncell;
                      
                      for (int j=0; j< Nmax; j++) {
                          int jloc = cell_list[s*Ncell*Ncell+aloc*Ncell + bloc][j];
                          //printf("%d %d\n",j,jloc);
                          if (jloc != -1 && jloc != i) {

                              if (type_data[jloc+s*N]==0) {
                          
                                double dr=0.0, dx[NDim];
                                for (int d=0; d<NDim;d++) {
                                    dx[d] = xyz_data[i+s*N][d] - xyz_data[jloc+s*N][d];
                                    apply_pbc(dx[d]);
                                    dr += dx[d]*dx[d];
                                }
                                
                                int bin = sqrt(dr)/rc*NBin;
                                if (bin < NBin) {
                                  G4[bin] += (input_tensor[i+s*N])*(input_tensor[jloc+s*N]);
                                  //G4[bin] += 1.0;
                                }

                              }

                          }
                      }
                  }
              }

              sum_pred += (input_tensor[i+s*N]);
              //sum_pred += 1.0;

            }
        }
    }

		double NsAvg_pred = sum_pred/((double)NS);

    printf("%f\n",NsAvg_pred);

	 // print
    pathPred = QString::fromStdString(folderOut);
    pathPred.append(QString("/G4_%1.dat").arg(QString::fromStdString(input_name)));

    QFile outfilePred4(pathPred);  
    outfilePred4.open(QIODevice::WriteOnly | QIODevice::Text );
    QTextStream outPred4(&outfilePred4);

		for (int bin=0; bin<NBin;bin++) {
      double scale_pred = boxL*boxL/(NS*3.1415926*((double) (DeltaBin*DeltaBin)*(2*bin+1)))/(NsAvg_pred*(NsAvg_pred-1.0));
      outPred4 <<  bin*DeltaBin << " " << G4[bin]*scale_pred << " " << G4[bin];
      outPred4 << "\n";
    }
    outfilePred4.close();
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