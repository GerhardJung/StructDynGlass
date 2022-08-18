from signal import SIGUSR1
import pandas as pd
import numpy as np
import os

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
from scipy.interpolate import LinearNDInterpolator
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

N=4096
NType=2
NArray = [3277, 819]
NameArray = ["type0", "type1"]
boxL=15.041

Temp=0.44
NT=10
NSTrain=30
NSTest_array=[20,20,20,20]
Tstart=7
Tstop=NT

dyn_array = ['MD']

# data to save particle predictions
"""NTsave=int(Tstop-Tstart)
Nmodes=len(dyn_array)
data_true = np.empty((NTsave*Nmodes,NSTest*N), dtype=float)
data_pred = np.empty((NTsave*Nmodes,NSTest*N), dtype=float)"""

# read time data
time_data = np.empty((NT+30), dtype=float)
iso_data = np.empty((NT+30), dtype=float)
f = open("times.dat", "r")
count = 1
for x in f:
	time_data[count]=x.split()[0]
	count+=1
	
	
"""axes1 = []
axes2 = []"""
def test_model(NLOC,Type,NSCount,NSTest):

	# read labels
	KA_labels_full = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T"+ str(Temp) +"_bapst/ml_labels_" + Type + ".csv",nrows=(NSTest)*NLOC+2, skiprows=[i for i in range(3,(NSCount+NSTrain)*NLOC+3)])
	KA_labels_full = KA_labels_full.iloc[: , :-1]
	
	#print(KA_labels_full)

	# read features
	KA_features_phys = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T" + str(Temp) +"_bapst/ml_struct_" + Type + ".csv",nrows=(NSTest)*NLOC+3, skiprows=[i for i in range(4,(NSCount+NSTrain)*NLOC+4)])
	
	#print(KA_features_phys)

	# extract positions
	XPOS = KA_features_phys.pop('NDim0')
	XPOS.drop([0,1,2], inplace=True)
	#print(XPOS)
	YPOS = KA_features_phys.pop('NDim1')
	YPOS.drop([0,1,2], inplace=True)
	ZPOS = KA_features_phys.pop('NDim2')
	ZPOS.drop([0,1,2], inplace=True)
	ID = KA_features_phys.pop('ID')
	ID.drop([0,1,2], inplace=True)
	KA_features_phys_length = KA_features_phys.values[0]
	KA_features_phys_mean = KA_features_phys.values[1]
	KA_features_phys_var = KA_features_phys.values[2]
	KA_features_phys.drop([0,1,2], inplace=True)
	KA_features_phys.index -= 3
	XPOS.index -= 3
	YPOS.index -= 3
	ZPOS.index -= 3
	ID.index -= 3
	#print(KA_features_phys)
	#print(ID)
	
	# rescale features
	labels=KA_features_phys.columns.values
	ind=0
	for var in KA_features_phys_var:
		#print(ind)
		KA_features_phys[labels[ind]] = KA_features_phys[labels[ind]].add(-KA_features_phys_mean[ind]).div(var)
		ind+=1
	length=KA_features_phys.shape[1]
	for ind in range(length-1,-1,-1):
		var=KA_features_phys_var[ind]
		#print(var)
		if var < 0.0001 :
			KA_features_phys.drop(KA_features_phys.columns[ind], axis=1, inplace=True)
			KA_features_phys_length=np.delete(KA_features_phys_length,ind)
		ind+=1


	count=0
	for d in dyn_array:
		sourceFile = open(PATH + "/pearson_{}_{}.dat".format(d, Type), 'w')
		for t in range(Tstart, Tstop):
			
			print ("Analyze ", Type, " for observable ", d, " at time ", t)

			# choose labels
			KA_labels = KA_labels_full["{}{}".format(d, t)]
			# extract rescale quantities
			Mean_True =KA_labels[0]	
			Var_True = KA_labels[1]	
			KA_labels = KA_labels[2:]	
			KA_labels.index -= 2
			#print(KA2D_labels)

			KA_features = KA_features_phys
			KA_labels_tanh = (tf.tanh(  (KA_labels-0.44)*20.0  ) + 1.0)/2.0
				
			length=KA_features.shape[1]
			#print(length)
			

			# load bapst data
			count=0
			y_out=np.zeros((NSTest)*NLOC)
			fp = open("./pred_propensity_bapst_{}.dat".format(t-1), 'r')
			for i, line in enumerate(fp):
				if i >= (NSCount+NSTrain)*N and i < (NSCount+NSTrain+NSTest)*N :
					arr=line.split()
					#print(arr[0])
					if arr[0] == "0" :
						y_out[count] = arr[1]
						count+=1
					
			# load ML model
			def loss_carrier():
				def loss(y_true, y_pred):
					loss = 1.0		
					return loss
				return loss
			KA_model = load_model("./models_mlp/mlpparams_{}{}_{}.dat".format(d,t, Type), custom_objects={'loss': loss_carrier()})
			   
			# evaluate model
			y_out_new=KA_model.predict_on_batch(KA_features)[:,0]	
			
				
					
			#print(y_out)
			#print(y_out_loc)
			#print(KA_labels)
			y_out_tanh=(tf.tanh(  (y_out-0.44)*20.0  ) + 1.0)/2.0
			y_out_tanh_new=(tf.tanh(  (y_out_new-0.44)*20.0  ) + 1.0)/2.0

			corr, _ = pearsonr(y_out, KA_labels)
			coef, _ = spearmanr(y_out, KA_labels)
			
			corr_new, _ = pearsonr(y_out_new, KA_labels)
			coef_new, _ = spearmanr(y_out_new, KA_labels)
			
			# save data
			"""for s in range(NSTest):
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
					count_loc += 1"""
			#print(data_pred[count])
			
			# calc chi_4
			var = mean = 0
			var_single = mean_single = 0
			for s in range(NSTest):
				qtot_true = 0
				labels = KA_labels_tanh[s*NLOC:(s+1)*NLOC].numpy()
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
			#print((var-mean*mean)/NLOC,var_single-mean_single*mean_single)
			chi4_true=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
			chi4_true_standard=(var-mean*mean)/NLOC
			var_true = (var_single-mean_single*mean_single)
			var = mean = 0
			var_single = mean_single = 0
			for s in range(NSTest):
				qtot_pred = 0
				labels = y_out_tanh[s*NLOC:(s+1)*NLOC].numpy()
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
			chi4_pred=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
			chi4_pred_standard=(var-mean*mean)/NLOC
			
			var = mean = 0
			var_single = mean_single = 0
			for s in range(NSTest):
				qtot_pred = 0
				labels = y_out_tanh_new[s*NLOC:(s+1)*NLOC].numpy()
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
			chi4_pred_new=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
			chi4_pred_standard_new=(var-mean*mean)/NLOC
			var_pred = (var_single-mean_single*mean_single)
			
			print (time_data[t], corr,coef,corr_new,coef_new, chi4_true, chi4_pred,chi4_pred_new, chi4_true_standard, chi4_pred_standard, chi4_pred_standard_new)
			#print (time_data[t],corr_new,coef_new, chi4_true, chi4_pred_new, chi4_true_standard,  chi4_pred_standard_new, var_true, var_pred)
			print (time_data[t], corr,coef,corr_new,coef_new, chi4_true, chi4_pred,chi4_pred_new, chi4_true_standard, chi4_pred_standard, chi4_pred_standard_new, file = sourceFile)
			
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
				y_outm2 = y_out[KA_labels<0.2]
				y_outm1 = y_out[(KA_labels<0.4) & (KA_labels>0.2)]
				y_out0 = y_out[(KA_labels<0.6) & (KA_labels>0.4)]
				y_outp1 = y_out[(KA_labels<0.8) & (KA_labels>0.6)]
				y_outp2 = y_out[(KA_labels<=1.0) & (KA_labels>0.8)]
				bins = np.linspace(0.75, 1, 100)
			
			#plt.hist([y_outrest,y_outlow,y_outhigh], bins,  alpha=0.5, label=['model rest','model low','model high'], stacked=True)
			plt.figure(1)
			plt.clf()
			plt.hist(KA_labels, bins,  alpha=0.5, label='true')
			plt.hist([y_outm2,y_outm1,y_out0,y_outp1,y_outp2], bins,  alpha=0.5, label=['model <0.2','model <0.4','model <0.6','model <0.8','model <1.0'], histtype='step', stacked=True, fill=False)
			plt.legend(loc='upper left')
			plt.title("{} t={}".format(d, time_data[t]))
			plt.xlabel("{}".format(d))
			plt.ylabel("PDF")
			plt.savefig(PATH + "/histogram_{}{}_{}.png".format(d, t, Type))
			
			plt.figure(2)
			plt.clf()
			plt.scatter(KA_labels, y_out,s=1,  c="g", alpha=0.5)
			plt.title("{} t={}".format(d, time_data[t]))
			plt.xlabel("true")
			plt.ylabel("pred")
			plt.savefig(PATH + "/scatter_{}{}_{}.png".format(d, t, Type))
			
			
			# S_4
			qmin = 2.0 * 3.1415926 / boxL
			Nq = 40
			
			# calc density profiles
			bVec_x = qmin*XPOS
			bVec_y = qmin*YPOS
			bVec_z = qmin*ZPOS
			#print(bVec)
			resVec = []
			for k in range (3) :
				if k == 0 :
					bVec = bVec_x
				if k == 1 :
					bVec = bVec_y
				if k == 2 :
					bVec = bVec_z
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
			
			redResVec_true = []
			redResVec_pred = []
			#print (y_out)
			#print(KA2D_labels_test)
			#print(resVec[0])
			
			for res in resVec :
				for s in range(NSTest):
					redResVec_true.append(np.sum(res[NLOC*(s):NLOC*((s)+1)]*(KA_labels_tanh[NLOC*(s):NLOC*((s)+1)])))
					redResVec_pred.append(np.sum(res[NLOC*(s):NLOC*((s)+1)]*(y_out_tanh[NLOC*(s):NLOC*((s)+1)])))
			
			# calc structure factor
			strucFac_true = np.zeros(Nq)
			strucFac_pred = np.zeros(Nq)
			for k in range (3) :
				for m in range (Nq) :
					for s in range(NSTest):
						c_true = redResVec_true[2*k*Nq*NSTest + 2*m*NSTest +s]	
						s_true = redResVec_true[2*k*Nq*NSTest + 2*m*NSTest + NSTest +s]
						strucFac_true[m] += c_true*c_true + s_true*s_true	
						c_pred = redResVec_pred[2*k*Nq*NSTest + 2*m*NSTest +s]	
						s_pred = redResVec_pred[2*k*Nq*NSTest + 2*m*NSTest + NSTest +s]	
						strucFac_pred[m] += c_pred*c_pred + s_pred*s_pred
					
			
			# print
			sourceFileS4 = open(PATH + "/S4_{}{}_{}.dat".format(d,t, Type), 'w')
			for m in range (Nq) :
				#print ((m+1)*qmin, strucFac_true[m]/(2.0*NLOC*NSTest), strucFac_pred[m]/(2.0*NLOC*NSTest))
				print ((m+1)*qmin, strucFac_true[m]/(3.0*NLOC*NSTest), strucFac_pred[m]/(3.0*NLOC*NSTest),file = sourceFileS4)
			
			sourceFileS4.close()
			
			# g_4
			"""NBin = 50
			DeltaBin = 0.2
			G4 = np.zeros(NBin)
			G4_pred = np.zeros(NBin)
			sum_true = 0
			sum_pred=0
			#G4_count = np.zeros(50)
			
			# calc distance matrices, norms and true label results	
			for s in range(1):
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
			size=30.0
			#iso_true=np.mean(KA2D_labels_test)
			#iso_pred=np.mean(y_out)
			vmin = 0.2
			vmax = 0.8
			scales = np.linspace(vmin, vmax, 10)
			cmap = plt.get_cmap('coolwarm_r')
			norm = plt.Normalize(scales.min(), scales.max())
			
			size_ml = size*y_out_tanh[-NLOC:]
			size_ml_new = size*y_out_tanh_new[-NLOC:]
			size_true = size*KA_labels_tanh[-NLOC:]
			
			sloc=4
			points = (XPOS[-NLOC-sloc*NLOC:-sloc*NLOC], YPOS[-NLOC-sloc*NLOC:-sloc*NLOC],ZPOS[-NLOC-sloc*NLOC:-sloc*NLOC])
			x = np.linspace(min(XPOS[-NLOC-sloc*NLOC:-sloc*NLOC])+0.1, max(XPOS[-NLOC-sloc*NLOC:-sloc*NLOC])-0.1, num=150)
			y = np.linspace(min(YPOS[-NLOC-sloc*NLOC:-sloc*NLOC])+0.1, max(YPOS[-NLOC-sloc*NLOC:-sloc*NLOC])-0.1, num=150)
			z = 0.0
			X, Y, Z = np.meshgrid(x, y,z)
			X2D, Y2D = np.meshgrid(x, y)
			
			vmin=0.2
			vmax=0.8
			
			plt.figure(3)
			plt.clf()
			linInterML= LinearNDInterpolator(points, y_out[-NLOC-sloc*NLOC:-sloc*NLOC])
			VALML = linInterML(X, Y, Z)[:-1, :-1,0]
			plt.pcolormesh(X2D, Y2D, VALML, cmap='coolwarm', vmin=vmin, vmax=vmax)
			plt.colorbar() # Color Bar
			plt.savefig(PATH + "/conf_{}{}_ml2D.png".format(d, t))
			
			plt.figure(4)
			plt.clf()
			linInterML= LinearNDInterpolator(points, y_out_new[-NLOC-sloc*NLOC:-sloc*NLOC])
			VALML = linInterML(X, Y, Z)[:-1, :-1,0]
			plt.pcolormesh(X2D, Y2D, VALML, cmap='coolwarm', vmin=vmin, vmax=vmax)
			plt.colorbar() # Color Bar
			plt.savefig(PATH + "/conf_{}{}_new2D.png".format(d, t))
			
			plt.figure(5)
			plt.clf()
			linInterML= LinearNDInterpolator(points, KA_labels[-NLOC-sloc*NLOC:-sloc*NLOC])
			VALML = linInterML(X, Y, Z)[:-1, :-1,0]
			plt.pcolormesh(X2D, Y2D, VALML, cmap='coolwarm', vmin=vmin, vmax=vmax)
			plt.colorbar() # Color Bar
			plt.savefig(PATH + "/conf_{}{}_md2D.png".format(d, t))	
			
			#print(VALML)
			
			shift = 6
			if Type == 'type0' :
				ind =  (t-Tstart) +shift
				fig = plt.figure(ind)
				ax = fig.add_subplot(111, projection='3d')
				
				cube = ax.scatter(XPOS[-NLOC:], YPOS[-NLOC:],ZPOS[-NLOC:], zdir='z', c=y_out[-NLOC:],cmap='coolwarm', s = size_ml, vmin=vmin, vmax=vmax)  # Plot the cube	
				ax.axes.xaxis.set_ticks([])
				ax.axes.yaxis.set_ticks([])
				ax.axes.zaxis.set_ticks([])
				plt.tight_layout()

				#cbar = fig.colorbar(cube, shrink=0.6, aspect=5) 
				
				fig.savefig(PATH + "/conf_{}{}_ml.png".format(d, t))
				
				ind =  (t-Tstart) +shift + NT
				fig = plt.figure(ind)
				ax = fig.add_subplot(111, projection='3d')

				
				cube = ax.scatter(XPOS[-NLOC:], YPOS[-NLOC:],ZPOS[-NLOC:], zdir='z', c=y_out_new[-NLOC:],cmap='coolwarm', s = size_ml_new, vmin=vmin, vmax=vmax)  # Plot the cube 	
				ax.axes.xaxis.set_ticks([])
				ax.axes.yaxis.set_ticks([])
				ax.axes.zaxis.set_ticks([])
				#cbar = fig.colorbar(cube, shrink=0.6, aspect=5) 
				
				fig.savefig(PATH + "/conf_{}{}_new_ml.png".format(d, t))
				
				ind =  (t-Tstart) +shift + 2*NT
				fig = plt.figure(ind)
				ax = fig.add_subplot(111, projection='3d')
				
				cube = ax.scatter(XPOS[-NLOC:], YPOS[-NLOC:],ZPOS[-NLOC:], zdir='z', c=KA_labels[-NLOC:],cmap='coolwarm', s = size_true, vmin=vmin, vmax=vmax)  # Plot the cube  
				ax.axes.xaxis.set_ticks([])
				ax.axes.yaxis.set_ticks([])
				ax.axes.zaxis.set_ticks([]) 
				cbar = fig.colorbar(cube, shrink=0.6, aspect=5) 	
				
				fig.savefig(PATH + "/conf_{}{}_true.png".format(d, t))			
				
				
			"""if Type == 'type1' :
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


		sourceFile.close()
					
# create folder for learned models
NSCount=0
for NSTest in NSTest_array:
	PATH='./results_mlp_' + str(NSCount)
	isExist = os.path.exists(PATH)
	if not isExist:
	  os.makedirs(PATH)
						
	# execute for the three different types					
	for type in range(1):
		test_model(NArray[type],NameArray[type],NSCount, NSTest)
		
	NSCount += NSTest
		

