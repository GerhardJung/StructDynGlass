from signal import SIGUSR1
import pandas as pd
import numpy as np
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# fitting
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge

# Make numpy values easier to read.
np.set_printoptions(precision=3, suppress=True)

# plotting:
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams['figure.dpi']=300 # highres display

from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy import sparse
import statistics

import tensorflow as tf
from tensorflow.python.ops import math_ops
from tensorflow.keras import layers
from keras.models import load_model
from keras.layers import Dropout
from tensorflow.keras import optimizers
from tensorflow.keras import initializers
from tensorflow.keras import regularizers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

N=1290
NType=3
NArray = [600, 330, 360]
NameArray = ["type0", "type1", "type2"]
boxL=32.8962

Temp=0.231
NT=41
NSTrain=200
NSTest_array=[25,25,25,25]
Tstart=40
Tstop=41

dyn_array = ['BB']

# data to save particle predictions
"""NTsave=int(Tstop-Tstart)
Nmodes=len(dyn_array)
data_true = np.empty((NTsave*Nmodes,NSTest*N), dtype=float)
data_pred = np.empty((NTsave*Nmodes,NSTest*N), dtype=float)"""

# read time data
time_data = np.empty((NT+30), dtype=float)
iso_data = np.empty((NT+30), dtype=float)
f = open("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"/isoconf_predictability_BB.dat", "r")
count = 1
for x in f:
	time_data[count]=x.split()[0]
	iso_data[count]=x.split()[1]
	count+=1
	
# create folder for learned models
PATH='./results_mlp'
isExist = os.path.exists(PATH)
if not isExist:
  os.makedirs(PATH)
	
axes1 = []
axes2 = []
def test_model(NLOC,Type,NSCount,NSTest):

	# read labels
	KA2D_labels_full = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T"+ str(Temp) +"/ml_labels_" + Type + ".csv",nrows=(NSTest)*NLOC+2, skiprows=[i for i in range(3,(NSCount+NSTrain)*NLOC+3)])
	KA2D_labels_full = KA2D_labels_full.iloc[: , :-1]

	# read features
	KA2D_features_phys = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"/ml_struct_" + Type + ".csv",nrows=(NSTest)*NLOC+3, skiprows=[i for i in range(4,(NSCount+NSTrain)*NLOC+4)])

	XPOS = KA2D_features_phys.pop('NDim0')
	XPOS.drop([0,1,2], inplace=True)
	YPOS = KA2D_features_phys.pop('NDim1')
	YPOS.drop([0,1,2], inplace=True)
	ID = KA2D_features_phys.pop('ID')
	ID.drop([0,1,2], inplace=True)
	KA2D_features_phys_length = KA2D_features_phys.values[0]
	KA2D_features_phys_mean = KA2D_features_phys.values[1]
	KA2D_features_phys_var = KA2D_features_phys.values[2]
	KA2D_features_phys.drop([0,1,2], inplace=True)
	KA2D_features_phys.index -= 3
	XPOS.index -= 3
	YPOS.index -= 3
	ID.index -= 3
	print(KA2D_features_phys)
	#print(ID)
	
	# rescale features
	labels=KA2D_features_phys.columns.values
	ind=0
	for var in KA2D_features_phys_var:
		#print(ind)
		KA2D_features_phys[labels[ind]] = KA2D_features_phys[labels[ind]].add(-KA2D_features_phys_mean[ind]).div(var)
		ind+=1
	length=KA2D_features_phys.shape[1]
	for ind in range(length-1,-1,-1):
		var=KA2D_features_phys_var[ind]
		#print(var)
		if var < 0.0001 :
			KA2D_features_phys.drop(KA2D_features_phys.columns[ind], axis=1, inplace=True)
			KA2D_features_phys_length=np.delete(KA2D_features_phys_length,ind)
		ind+=1
	#length=KA2D_features_phys.shape[1]
	#for ind in range(length-1,-1,-1):
	#	lengthscale=KA2D_features_phys_length[ind]
	#	if lengthscale > 5.75 :
	#		KA2D_features_phys.drop(KA2D_features_phys.columns[ind], axis=1, inplace=True)
	#		KA2D_features_phys_length=np.delete(KA2D_features_phys_length,ind)
	#	ind+=1

	count=0
	for d in dyn_array:
	
		Nmod = 1
		for mod in range(Nmod):
			#sourceFile = open(PATH + "/pearson_{}_{}_{}.dat".format(d, Type,mod), 'w')
																											
			for t in range(Tstart, Tstop):
				
				print ("Analyze ", Type, " for observable ", d, " at time ", t)

				# choose labels
				KA2D_labels_test = KA2D_labels_full["{}{}".format(d, t)]
				KA2D_features_test = KA2D_features_phys
				# extract rescale quantities
				Mean_True =KA2D_labels_test[0]	
				Var_True = KA2D_labels_test[1]	
				KA2D_labels_test = KA2D_labels_test[2:]	
				KA2D_labels_test.index -= 2
				if d == "BB" :
					KA2D_labels_test_mod = 1.0 - KA2D_labels_test
				if d == "MD" :
					KA2D_labels_test_mod = (tf.tanh(  (KA2D_labels_test-0.4)*20.0  ) + 1.0)/2.0
				#print(KA2D_labels)
				
				print(KA2D_features_test)
				print(KA2D_labels_test)
				print(KA2D_labels_test_mod)
					
				length=KA2D_features_test.shape[1]
				#print(length)
				
				# shuffle arrays
				#KA2D_features = KA2D_features.sample(frac=1,random_state=12345)
				#KA2D_labels = KA2D_labels.sample(frac=1,random_state=12345)
				#select labels
				#NumSelect = int((NLOC*NSTrain)/5)
				#smallest_index = KA2D_labels.nsmallest(NumSelect).index
				#largest_index = KA2D_labels.nlargest(NumSelect).index
				#select_index=smallest_index.union(largest_index)
				#print (largest_index)
				#KA2D_labels = KA2D_labels.iloc[select_index-2]
				#print(KA2D_labels)
				#KA2D_features = KA2D_features.iloc[select_index-2]
				
				# load ML model
				def loss_carrier():
					def loss(y_true, y_pred):
						loss = 1.0		
						return loss
					return loss
				KA2D_model = load_model("./models_mlp/mlpparams_{}{}_{}_{}.dat".format(d,t, Type,mod), custom_objects={'loss': loss_carrier()})
				   
				# evaluate model
				y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				if d == "BB" :
					y_out_mod = 1.0 - y_out
				if d == "MD" :
					y_out_mod = (tf.tanh(  (y_out-0.4)*20.0  ) + 1.0)/2.0
				#y_out *= Var_True
				#y_out += Mean_True
				#y_out=np.where(y_out > 1.0, 1.0, y_out)
				print(KA2D_labels_test.shape,y_out.shape)
				corr, _ = pearsonr(y_out, KA2D_labels_test)
				coef, _ = spearmanr(y_out, KA2D_labels_test)
				
				# save data
				for s in range(NSTest):
					KA2D_labels_test_loc = KA2D_labels_test[NLOC*s:NLOC*(s+1)].reset_index(drop=True)
					y_out_loc = y_out[NLOC*s:NLOC*(s+1)]
					ID_loc = ID[NLOC*(s+NSTrain):NLOC*(s+NSTrain+1)].astype(int)
					#if s == 0 :
					#	print (ID_loc)
					count_loc = 0
					for id in ID_loc :
						#if s == 0 :
						#	print (id)
						data_pred[count,id+s*N] = y_out_loc[count_loc]
						data_true[count,id+s*N] = KA2D_labels_test_loc[count_loc]
						count_loc += 1
				#print(data_pred[count])
				
				# calc chi_4
				var = mean = 0
				var_single = mean_single = 0
				for s in range(NSTest):
					qtot_true = 0
					labels = KA2D_labels_test_mod[s*NLOC:(s+1)*NLOC]
					for x in labels :
						qtot_true+=x
						mean_single += x
						var_single += x*x
					#print(s,qtot_pred)
					mean += qtot_true
					var += qtot_true*qtot_true
				mean/=NSTest
				var/=NSTest
				mean_single/=NSTest*NLOC
				var_single/=NSTest*NLOC
				mean_single = var_single
				var_single = var_single
				#print((var-mean*mean)/NLOC,var_single-mean_single*mean_single)
				chi4_true=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
				chi4_true_standard=(var-mean*mean)/NLOC
				var = mean = 0
				var_single_ml = mean_single_ml = 0
				for s in range(NSTest):
					qtot_pred = 0
					labels = y_out_mod[s*NLOC:(s+1)*NLOC]
					for x in labels :
						qtot_pred+=x
						mean_single_ml += x
						var_single_ml += x*x
					#print(s,qtot_pred)
					mean += qtot_pred
					var += qtot_pred*qtot_pred
				mean/=NSTest
				var/=NSTest
				mean_single_ml/=NSTest*NLOC
				var_single_ml/=NSTest*NLOC
				mean_single_ml = var_single_ml
				var_single_ml = var_single_ml
				chi4_pred=(var-mean*mean)/NLOC/(var_single_ml-mean_single_ml*mean_single_ml)
				chi4_pred_standard=(var-mean*mean)/NLOC
				
				
				print (time_data[t], corr,coef, chi4_true, chi4_pred, chi4_true_standard, chi4_pred_standard, var_single-mean_single*mean_single, var_single_ml-mean_single_ml*mean_single_ml)
				#print (time_data[t], corr,coef, chi4_true, chi4_pred, chi4_true_standard, chi4_pred_standard, file = sourceFile)
				
				# print histograms
				if d=='LOG(FRES)' :
					y_outm2 = y_out[KA2D_labels_test<-2.0]
					y_outm1 = y_out[(KA2D_labels_test<-1.0) & (KA2D_labels_test>-2.0)]
					y_out0 = y_out[(KA2D_labels_test<0.0) & (KA2D_labels_test>-1.0)]
					y_outp1 = y_out[(KA2D_labels_test<1.0) & (KA2D_labels_test>0)]
					y_outp2 = y_out[(KA2D_labels_test<2.0) & (KA2D_labels_test>1.0)]
					y_outp3 =y_out[(KA2D_labels_test>2.0)]
					bins = np.linspace(-4, 4, 100)
				if d=='MD' or d=='BB':
					y_outm2 = y_out[KA2D_labels_test<0.2]
					y_outm1 = y_out[(KA2D_labels_test<0.4) & (KA2D_labels_test>0.2)]
					y_out0 = y_out[(KA2D_labels_test<0.6) & (KA2D_labels_test>0.4)]
					y_outp1 = y_out[(KA2D_labels_test<0.8) & (KA2D_labels_test>0.6)]
					y_outp2 = y_out[(KA2D_labels_test<=1.0) & (KA2D_labels_test>0.8)]
					bins = np.linspace(0.0, 1, 100)
				
				#plt.hist([y_outrest,y_outlow,y_outhigh], bins,  alpha=0.5, label=['model rest','model low','model high'], stacked=True)
				plt.figure(1)
				plt.clf()
				hist_true=plt.hist(KA2D_labels_test, bins,  alpha=0.5, label='true')
				hist_pred=plt.hist([y_outm2,y_outm1,y_out0,y_outp1,y_outp2], bins,  alpha=0.5, label=['model <0.2','model <0.4','model <0.6','model <0.8','model <1.0'], histtype='step', stacked=True, fill=False)
				pd.DataFrame({'x_upper':hist_true[1][1:], 'y': hist_true[0]}).to_csv(PATH + "/histogram_datatrue_{}{}_{}_{}.dat".format(d, t, Type, mod),sep=' ')
				pd.DataFrame({'x_upper':hist_pred[1][1:], 'y02': hist_pred[0][0], 'y04': hist_pred[0][1], 'y06': hist_pred[0][2], 'y08': hist_pred[0][3], 'y': hist_pred[0][4]}).to_csv(PATH + "/histogram_datapred_{}{}_{}_{}.dat".format(d, t, Type, mod),sep=' ')
				plt.legend(loc='upper left')
				plt.title("{} t={}".format(d, time_data[t]))
				plt.xlabel("{}".format(d))
				plt.ylabel("PDF")
				plt.savefig(PATH + "/histogram_{}{}_{}.png".format(d, t, Type))
				
				plt.figure(2)
				plt.clf()
				plt.scatter(KA2D_labels_test, y_out,s=1,  c="g", alpha=0.5)
				plt.title("{} t={}".format(d, time_data[t]))
				plt.xlabel("true")
				plt.ylabel("pred")
				plt.savefig(PATH + "/scatter_{}{}_{}.png".format(d, t, Type))
				
				
				# S_4
				qmin = 2.0 * 3.1415926 / boxL
				Nq = 10
				
				# calc density profiles
				bVec_x = qmin*XPOS[:]
				bVec_y = qmin*YPOS[:]
				#print(bVec)
				resVec = []
				for k in range (2) :
					if k == 0 :
						bVec = bVec_x
					if k == 1 :
						bVec = bVec_y
					for m in range (Nq) :
						if m == 0:
							c = np.cos(bVec)
							s = np.sin(bVec)
							c0 = c
						elif m == 1 :
							c1 = c
							s1 = s
							c = 2.0 * c0*c1 - 1.0
							s = 2.0 * c0*s1
						else :
							c2 = c1
							s2 = s1
							c1 = c
							s1 = s
							c = 2.0 * c0 * c1 - c2
							s = 2.0 * c0 * s1 - s2
						resVec_loc_c = c
						resVec_loc_s = s
						resVec.append(resVec_loc_c)
						resVec.append(resVec_loc_s)
				
				Re1_true = []
				Re2_true = []
				Im1_true = []
				Im2_true = []
				Re1_pred = []
				Re2_pred = []
				Im1_pred = []
				Im2_pred = []
				
				for Nx in range (Nq+1) :
					for Ny in range (Nq+1) :
						for s in range(NSTest):
						
							#if Nx == 0 or Ny == 0 :
								if Nx == 0 :
									Rex = np.ones(NLOC)
									Imx = np.zeros(NLOC)
								else :
									Rex = resVec[ 2*(Nx-1)][NLOC*(s):NLOC*((s)+1)]
									Imx = resVec[ 2*(Nx-1) + 1][NLOC*(s):NLOC*((s)+1)]
								if Ny == 0 :	
									Rey = np.ones(NLOC)
									Imy = np.zeros(NLOC)
								else :
									Rey = resVec[2*Nq + 2*(Ny-1)][NLOC*(s):NLOC*((s)+1)]
									Imy = resVec[2*Nq + 2*(Ny-1) + 1][NLOC*(s):NLOC*((s)+1)]	

								Re1_true.append(np.sum(Rex*Rey*KA2D_labels_test_mod[NLOC*(s):NLOC*((s)+1)]))
								Re2_true.append(np.sum(-Imx*Imy*KA2D_labels_test_mod[NLOC*(s):NLOC*((s)+1)]))
								Im1_true.append(np.sum(Rex*Imy*KA2D_labels_test_mod[NLOC*(s):NLOC*((s)+1)]))
								Im2_true.append(np.sum(Imx*Rey*KA2D_labels_test_mod[NLOC*(s):NLOC*((s)+1)]))
								
								Re1_pred.append(np.sum(Rex*Rey*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								Re2_pred.append(np.sum(-Imx*Imy*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								Im1_pred.append(np.sum(Rex*Imy*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								Im2_pred.append(np.sum(Imx*Rey*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
							#else :
								#Re1_true.append(0.0)
								#Re2_true.append(0.0)
								#Im1_true.append(0.0)
								#Im2_true.append(0.0)
								
								#Re1_pred.append(0.0)
								#Re2_pred.append(0.0)
								#Im1_pred.append(0.0)
								#Im2_pred.append(0.0)
				
				# calc structure factor
				strucFac_true = np.zeros( (Nq+1)*(Nq+1))
				strucFac_pred = np.zeros((Nq+1)*(Nq+1))
				for Nx in range (Nq+1) :
					for Ny in range (Nq+1) :
						for s in range(NSTest):
							c_true = Re1_true[ (Nx*(Nq+1)+Ny)*NSTest + s] + Re2_true[ (Nx*(Nq+1)+Ny)*NSTest + s] 
							s_true = Im1_true[ (Nx*(Nq+1)+Ny)*NSTest + s] + Im2_true[ (Nx*(Nq+1)+Ny)*NSTest + s] 
							strucFac_true[Nx*(Nq+1)+Ny] += c_true*c_true + s_true*s_true	
								
							c_pred = Re1_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] + Re2_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] 
							s_pred = Im1_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] + Im2_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] 
							strucFac_pred[Nx*(Nq+1)+Ny] += c_pred*c_pred + s_pred*s_pred
						
				
				# print
				sourceFileS4 = open(PATH + "/S4_{}{}_{}_{}.dat".format(d,t, Type, mod), 'w')
				for Nx in range (Nq+1) :
					for Ny in range (Nq+1) :
						if Nx != 0 or Ny != 0:
							qx = (Nx)*qmin
							qy = (Ny)*qmin
							print (   np.sqrt(qx*qx + qy*qy) , strucFac_true[Nx*(Nq+1)+Ny]/(NLOC*NSTest), strucFac_pred[Nx*(Nq+1)+Ny]/(NLOC*NSTest), strucFac_true[Nx*(Nq+1)+Ny]/(NLOC*NSTest)/(var_single-mean_single*mean_single), strucFac_pred[Nx*(Nq+1)+Ny]/(NLOC*NSTest)/(var_single_ml-mean_single_ml*mean_single_ml),file = sourceFileS4)
				
				sourceFileS4.close()
				
				# g_4
				NBin = 200
				DeltaBin = 0.05
				G4 = np.zeros(NBin)
				G4_pred = np.zeros(NBin)
				sum_true = 0
				sum_pred=0
				#G4_count = np.zeros(50)
				
				# calc distance matrices, norms and true label results	
				for s in range(NSTest):
					#print(s)
					# calc distances
					positions = np.stack((XPOS[NLOC*(s):NLOC*((s)+1)],YPOS[NLOC*(s):NLOC*((s)+1)]),axis=1 )
					dist_nd_sq = np.zeros(NLOC * (NLOC - 1) // 2)  # to match the result of pdist
					for dim in range(2):
						pos_1d = positions[:, dim][:, np.newaxis]  # shape (N, 1)
						dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
						dist_1d[dist_1d > boxL * 0.5] -= boxL
						dist_1d[dist_1d < -boxL * 0.5] += boxL
						#print (dist_nd_sq,dist_1d)
						dist_nd_sq += dist_1d ** 2  # d^2 = dx^2 + dy^2 + dz^2
					#dist_nd = squareform(np.sqrt(dist_nd_sq))
					dist_nd = np.sqrt(dist_nd_sq)/DeltaBin
					bin_nd = dist_nd.astype(int)
					bin_nd[bin_nd >= NBin] = 0
					bin_nd = squareform (bin_nd)
					bin_nd = sparse.csr_matrix(bin_nd)
					
					KA2D_labels_test_loc = 1.0-KA2D_labels_test[NLOC*s:NLOC*(s+1)].reset_index(drop=True)
					#KA2D_labels_test_loc=np.where(KA2D_labels_test_loc > 0.2, 1.0, KA2D_labels_test_loc)
					#KA2D_labels_test_loc=np.where(KA2D_labels_test_loc <= 0.2, 0.0, KA2D_labels_test_loc)
					y_out_loc =1.0-y_out[NLOC*s:NLOC*(s+1)]
					#y_out_loc=np.where(y_out_loc > 0.2, 1.0, y_out_loc)
					#y_out_loc=np.where(y_out_loc <= 0.2, 0.0, y_out_loc)
					sum_true += np.sum(KA2D_labels_test_loc)
					sum_pred += np.sum(y_out_loc)
					#print(np.sum(KA2D_labels_test_loc),np.sum(y_out_loc))
					#print(KA2D_labels_test_loc)
					#print(y_out_loc)
					#print(KA2D_labels_test_loc)
					for i, j in zip(*bin_nd.nonzero()):
						G4[bin_nd[i,j]] += KA2D_labels_test_loc[i]*KA2D_labels_test_loc[j]
						G4_pred[bin_nd[i,j]] += y_out_loc[i]*y_out_loc[j]
						#G4[bin_nd[i,j]] += 1.0
						#G4_count[bin_nd[i,j]] += 1.0
						#print(row,col,bin_nd[row, col])
						
					#print (bin_nd)

			
				NsAvg_true = sum_true/NSTest
				NsAvg_pred = sum_pred/NSTest
				#NsAvg = NLOC
				sourceFileG4 = open(PATH + "/G4_{}{}_{}_{}.dat".format(d,t, Type, mod), 'w')
				for bin in range(0,NBin) :
					scale = boxL*boxL/(NSTest*3.1415926*DeltaBin*DeltaBin*(2*bin+1))/(NsAvg_true*(NsAvg_true))
					scale_pred = boxL*boxL/(NSTest*3.1415926*DeltaBin*DeltaBin*(2*bin+1))/(NsAvg_pred*(NsAvg_pred))
					#scale = scale_pred = boxL*boxL/(3.1415926*DeltaBin*(2*bin+1))
					#print (bin*DeltaBin, G4[bin]*scale, G4_pred[bin]*scale_pred,G4[bin],G4_pred[bin])
					print (bin*DeltaBin, G4[bin]*scale, G4_pred[bin]*scale_pred,G4[bin],G4_pred[bin],file = sourceFileG4)
				sourceFileG4.close()
				
				# print configurations
				"""size=30.0
				iso_true=np.mean(KA2D_labels_test)
				iso_pred=np.mean(y_out)
				vmin = 0.0
				vmax = 1.0
				scales = np.linspace(vmin, vmax, 10)
				cmap = plt.get_cmap('coolwarm_r')
				norm = plt.Normalize(scales.min(), scales.max())
				
				shift = 3
				if Type == 'type0' :
					ind =  (t-Tstart) +shift
					fig = plt.figure(ind, figsize=(10, 5))
					
					plt.title("{} t={}\n iso_pred={:.2f} iso_true={:.2f}".format(d, time_data[t], iso_pred, iso_true))
					plt.axis('off')
					ax1 = fig.add_subplot(1, 2, 1, aspect = "equal")
					ax1.set_xlabel("X")
					ax1.set_xlabel("Y")
					im1=ax1.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					axes1.append(ax1)
					
					ax2 = fig.add_subplot(1, 2, 2, aspect = "equal")
					ax2.set_xlabel("X")
					ax2.set_xlabel("Y")
					im2=ax2.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=KA2D_labels_test[-NLOC:],s=size, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					axes2.append(ax2)
					
					#divider1 = make_axes_locatable(ax1)
					#cax1 = divider1.append_axes("right", size="5%", pad=0.05)

					divider2 = make_axes_locatable(ax2)
					cax2 = divider2.append_axes("right", size="5%", pad=0.05)

					#Create and remove the colorbar for the first subplot
					#cbar1 = fig.colorbar(im1, cax = cax1)
					#fig.delaxes(fig.axes[1])

					#Create second colorbar
					cbar2 = fig.colorbar(im2, cax = cax2)

					#plt.tight_layout()					
					
					
				if Type == 'type1' :
					ind =  (t-Tstart) +shift
					print(ind)
					fig = plt.figure(ind)
					#ax1 = fig.add_subplot(1, 2, 1, aspect = "equal")
					im1=axes1[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size*0.88, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					
					#ax2 = fig.add_subplot(1, 2, 2, aspect = "equal")
					im2=axes2[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=KA2D_labels_test[-NLOC:],s=size*0.88, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
				
				if Type == 'type2' :
					ind =  (t-Tstart) +shift
					fig = plt.figure(ind)
					#ax1 = fig.add_subplot(1, 2, 1, aspect = "equal")
					im1=axes1[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size*0.94, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					
					#ax2 = fig.add_subplot(1, 2, 2, aspect = "equal")
					im2=axes2[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=KA2D_labels_test[-NLOC:],s=size*0.94, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					
					fig.savefig(PATH + "/conf_{}{}.png".format(d, t))
					
					count += 1"""


			#sourceFile.close()
					
# execute for the three different types
NSCount=0
for NSTest in NSTest_array:					
	for type in range(1):
	
		PATH='./results_mlp_1290_' + str(NSCount)
		isExist = os.path.exists(PATH)
		if not isExist:
		  os.makedirs(PATH)
							
		# execute for the three different types					
		for type in range(1):
			test_model(NArray[type],NameArray[type],NSCount, NSTest)
			
		NSCount += NSTest
	
			
# print predicted configurations
"""config_file = open(PATH +"/pred_propensity.dat", "w")
# header
config_file.write("ID ")
for d in dyn_array:
	for t in range(Tstart,Tstop):
		config_file.write("{}{}_PRED ".format(d,t))
		config_file.write("{}{}_TRUE ".format(d,t))
config_file.write("\n") 

for p in range(NSTest*N):
	config_file.write("%d " % p + NSTrain*N)
	count=0		
	for d in dyn_array:
		for t in range(Tstart,Tstop):
			config_file.write("%f " % data_pred[count,p])
			config_file.write("%f " % data_true[count,p])
			count += 1
	config_file.write("\n")"""

