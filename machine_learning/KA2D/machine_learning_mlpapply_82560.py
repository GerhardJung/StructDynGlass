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

N=82560
NType=3
NArray = [38400, 21120, 23040]
NameArray = ["type0", "type1", "type2"]
boxL=263.17

Temp=0.231
NT=41
NSTest_array=[20]
Tstart=27
Tstop=NT
INDEXarray=[1]

dyn_array = ['BB']



# data to save particle predictions
NPrint=20
NTsave=int(Tstop-Tstart)
Nmodes=len(dyn_array)
data_pred = np.empty((NTsave*Nmodes,NPrint*N), dtype=float)

# read time data

time_data = np.empty((NT+30), dtype=float)
iso_data = np.empty((NT+30), dtype=float)
f = open("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"/isoconf_predictability_BB.dat", "r")
count = 1
for x in f:
	time_data[count]=x.split()[0]
	iso_data[count]=x.split()[1]
	count+=1
	
"""axes1 = []
axes2 = []"""
def test_model(NLOC,Type,NSCount,NSTest,INDEX):

	# read features
	print("Start read data!")
	
	print ((NSCount+NSTest)*NLOC+3, NSCount*NLOC+3)
	
	KA2D_features_phys = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"_82560_"+str(INDEX)+"/ml_struct_" + Type + ".csv",nrows=(NSTest)*NLOC+3, skiprows=[i for i in range(4,NSCount*NLOC+4)])

	print("Finished read data!")

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
	print(XPOS)
	print(YPOS)
	print(KA2D_features_phys)
	
	print("Test0")
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
			sourceFile = open(PATH + "/pearson_{}_{}_{}.dat".format(d, Type,mod), 'w')
			for t in range(Tstart, Tstop):
				
				print ("Analyze ", Type, " for observable ", d, " at time ", t)


				KA2D_features_test = KA2D_features_phys
					
				length=KA2D_features_test.shape[1]
				#print(length)
				
				# load ML model
				def loss_carrier():
					def loss(y_true, y_pred):
						loss = 1.0		
						return loss
					return loss
					

				KA2D_model = load_model("./models_mlp/mlpparams_{}{}_{}_{}.dat".format(d,t, Type, mod), custom_objects={'loss': loss_carrier()})
				
				print("Test1")
				   
				# evaluate model
				y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				if d == "BB" :
					y_out_mod = 1.0 - y_out
				if d == "MD" :
					y_out_mod = (tf.tanh(  (y_out-0.4)*20.0  ) + 1.0)/2.0
				
				# save data
				if mod == 0 :
					for s in range(NSTest):
						y_out_loc = y_out_mod[NLOC*s:NLOC*(s+1)]
						ID_loc = ID[NLOC*(s):NLOC*(s+1)].astype(int)
						if s == 0 :
							print (ID_loc)
						count_loc = 0
						for id in ID_loc :
							#if s == 0 :
							#	print (id)
							data_pred[count,id+s*N+NSCount*N] = y_out_loc[count_loc]
							count_loc += 1
				#print(data_pred[count])
				
				print("Test2")
				
				# calc chi4
				var = mean = 0
				var_single = mean_single = 0
				for s in range(NSTest):
					qtot_pred = 0
					labels = y_out_mod[s*NLOC:(s+1)*NLOC]
					for x in labels :
						qtot_pred+=x
						mean_single += x
						var_single += x*x
					#print(s,qtot_pred)
					mean += qtot_pred
					var += qtot_pred*qtot_pred
				mean/=NSTest
				var/=NSTest
				mean_single/=NSTest*NLOC
				var_single/=NSTest*NLOC
				mean_single=mean_single
				var_single=var_single
				chi4_pred=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
				chi4_pred_standard = (var-mean*mean)/NLOC
				
				print (time_data[t], chi4_pred, chi4_pred_standard)
				print (time_data[t], chi4_pred, chi4_pred_standard, file = sourceFile)
				
				print("Var")
				print(np.var(y_out))	
					
				print("Mean")
				print(np.mean(y_out))
				
				"""# print histograms
				bins = np.linspace(0.75, 1, 100)
				plt.figure(1)
				plt.clf()
				plt.hist([y_out], bins,  alpha=0.5, histtype='step', stacked=True, fill=False)
				plt.legend(loc='upper left')
				plt.title("{} t={}".format(d, time_data[t]))
				plt.xlabel("{}".format(d))
				plt.ylabel("PDF")
				plt.savefig(PATH + "/histogram_{}{}_{}.png".format(d, t, Type))"""
				
				# S_4
				qmin = 2.0 * 3.1415926 / boxL
				Nq = 15
				
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
							#print (bVec)
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
				
				print (resVec[10])
				
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

								
								Re1_pred.append(np.sum(Rex*Rey*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								Re2_pred.append(np.sum(-Imx*Imy*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								Im1_pred.append(np.sum(Rex*Imy*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								Im2_pred.append(np.sum(Imx*Rey*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
								
								if s==0 and Nx==1 and Ny==2 :
									#print (Rex,Rey,y_out_mod[NLOC*(s):NLOC*((s)+1)])
									print(np.sum(Rex*Rey*y_out_mod[NLOC*(s):NLOC*((s)+1)]))
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
				strucFac_pred = np.zeros((Nq+1)*(Nq+1))
				for Nx in range (Nq+1) :
					for Ny in range (Nq+1) :
						for s in range(NSTest):
								
							c_pred = Re1_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] + Re2_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] 
							s_pred = Im1_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] + Im2_pred[ (Nx*(Nq+1)+Ny)*NSTest + s] 
							strucFac_pred[Nx*(Nq+1)+Ny] += c_pred*c_pred + s_pred*s_pred
							
						if Nx==1 and Ny==2 :
							print(strucFac_pred[Nx*(Nq+1)+Ny])
						
				
				# print
				sourceFileS4 = open(PATH + "/S4_{}{}_{}_{}.dat".format(d,t, Type, mod), 'w')
				for Nx in range (Nq+1) :
					for Ny in range (Nq+1) :
						if Nx != 0 or Ny != 0:
							qx = (Nx)*qmin
							qy = (Ny)*qmin
							print (   np.sqrt(qx*qx + qy*qy) , strucFac_pred[Nx*(Nq+1)+Ny]/(NLOC*NSTest), strucFac_pred[Nx*(Nq+1)+Ny]/(NLOC*NSTest)/(var_single-mean_single*mean_single),file = sourceFileS4)
				
				sourceFileS4.close()
				
				"""# g_4
				NBin = 50
				DeltaBin = 0.2
				G4 = np.zeros(NBin)
				G4_pred = np.zeros(NBin)
				sum_true = 0
				sum_pred=0
				#G4_count = np.zeros(50)
				
				# calc distance matrices, norms and true label results"""	
				"""for s in range(1):
					#print(s)
					# calc distances
					positions = np.stack((XPOS[NLOC*(s+NSTrain):NLOC*((s+NSTrain)+1)],YPOS[NLOC*(s+NSTrain):NLOC*((s+NSTrain)+1)]),axis=1 )
					dist_nd_sq = np.zeros(NLOC * (NLOC - 1) // 2)  # to match the result of pdist
					for dim in range(2):
						pos_1d = positions[:, dim][:, np.newaxis]  # shape (N, 1)
						dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
						dist_1d[dist_1d > boxL * 0.5] -= boxL
						dist_1d[dist_1d < -boxL * 0.5] += boxL
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
				sourceFileG4 = open(PATH + "/G4_{}{}_{}.dat".format(d,t, Type), 'w')
				for bin in range(0,NBin) :
					scale = boxL*boxL/(NSTest*3.1415926*DeltaBin*DeltaBin*(2*bin+1))/(NsAvg_true*(NsAvg_true))
					scale_pred = boxL*boxL/(NSTest*3.1415926*DeltaBin*DeltaBin*(2*bin+1))/(NsAvg_pred*(NsAvg_pred))
					#scale = scale_pred = boxL*boxL/(3.1415926*DeltaBin*(2*bin+1))
					#print (bin*DeltaBin, G4[bin]*scale, G4_pred[bin]*scale_pred,G4[bin],G4_pred[bin])
					print (bin*DeltaBin, G4[bin]*scale, G4_pred[bin]*scale_pred,G4[bin],G4_pred[bin],file = sourceFileG4)
				sourceFileG4.close()"""
				
				# print configurations
				"""size=30.0
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
					
					plt.title("{} t={}\n iso_pred={:.2f}".format(d, time_data[t], iso_pred))
					plt.axis('off')
					ax1 = fig.add_subplot(1, 2, 1, aspect = "equal")
					ax1.set_xlabel("X")
					ax1.set_xlabel("Y")
					im1=ax1.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					axes1.append(ax1)
					
					ax2 = fig.add_subplot(1, 2, 2, aspect = "equal")
					ax2.set_xlabel("X")
					ax2.set_xlabel("Y")
					im2=ax2.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
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
					im2=axes2[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size*0.88, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
				
				if Type == 'type2' :
					ind =  (t-Tstart) +shift
					fig = plt.figure(ind)
					#ax1 = fig.add_subplot(1, 2, 1, aspect = "equal")
					im1=axes1[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size*0.94, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					
					#ax2 = fig.add_subplot(1, 2, 2, aspect = "equal")
					im2=axes2[ind-shift].scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size*0.94, cmap='coolwarm_r', vmin=vmin, vmax=vmax)
					
					fig.savefig(PATH + "/conf_{}{}.png".format(d, t))"""
					
				count += 1


			sourceFile.close()
	
# create folder for learned models
for INDEX in INDEXarray:
	NSCount=0
	for NSTest in NSTest_array:
		print("INDEX", INDEX, "NSCount", NSCount)
		PATH='./results_mlp_82560_' + str(NSCount+(INDEX-1)*100)
		isExist = os.path.exists(PATH)
		if not isExist:
		  os.makedirs(PATH)
							
		# execute for the three different types					
		for type in range(3):
			test_model(NArray[type],NameArray[type],NSCount, NSTest, INDEX)
			
		NSCount += NSTest
	
			
# print predicted configurations
config_file = open("./eval_mlp/pred_propensity.dat", "w")
# header
config_file.write("ID ")
for d in dyn_array:
	for t in range(Tstart,Tstop):
		config_file.write("BB{}_PRED ".format(t))
config_file.write("\n") 

for p in range(NPrint*N):
	config_file.write("%d " % p )
	count=0		
	for d in dyn_array:
		for t in range(Tstart,Tstop):
			config_file.write("%f " % data_pred[count,p])
			count += 1
	config_file.write("\n")
		

