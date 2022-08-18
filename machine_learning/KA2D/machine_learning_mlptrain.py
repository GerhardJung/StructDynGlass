from signal import SIGUSR1
import pandas as pd
import numpy as np
from numpy import sqrt
from numpy import exp
import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# fitting
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge

# Make numpy values easier to read.
np.set_printoptions(precision=3, suppress=True)

# plotting:
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.dpi']=300 # highres display

from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import statistics

import tensorflow as tf
from tensorflow.python.ops import math_ops
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
from keras.models import load_model
from tensorflow.keras import layers
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
NSTrain=100
NSTest=100
Tstart=27
Tstop=NT

dyn_array = ['BB']

def train_model(NLOC,Type):

	# read simulation data
	# labels
	KA2D_labels_full = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T"+ str(Temp) +"/ml_labels_" + Type + ".csv",nrows=(NSTrain+NSTest)*NLOC+2)
	KA2D_labels_full = KA2D_labels_full.iloc[: , :-1]
	KA2D_labels_full_mean = KA2D_labels_full.values[0]
	KA2D_labels_full_var = KA2D_labels_full.values[1]
	KA2D_labels_full.drop([0,1], inplace=True)
	
	# rescale labels
	#labels=KA2D_labels_full.columns.values
	#ind=0
	#for var in KA2D_labels_full_var:
	#	KA2D_labels_full[labels[ind]] = KA2D_labels_full[labels[ind]].add(-KA2D_labels_full_mean[ind]).div(var)
	#	ind+=1

	# features
	KA2D_features_phys = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"/ml_struct_" + Type + ".csv",nrows=(NSTrain+NSTest)*NLOC+3)
	KA2D_features_phys = KA2D_features_phys.iloc[: , :-1]

	# extract positions
	XPOS = KA2D_features_phys.pop('NDim0')
	XPOS.drop([0,1,2], inplace=True)
	YPOS = KA2D_features_phys.pop('NDim1')
	YPOS.drop([0,1,2], inplace=True)
	KA2D_features_phys_length = KA2D_features_phys.values[0]
	KA2D_features_phys_mean = KA2D_features_phys.values[1]
	KA2D_features_phys_var = KA2D_features_phys.values[2]
	KA2D_features_phys.drop([0,1,2], inplace=True)
	
	# rescale features
	labels=KA2D_features_phys.columns.values
	ind=0
	for var in KA2D_features_phys_var:
		#print(ind)
		KA2D_features_phys[labels[ind]] = KA2D_features_phys[labels[ind]].add(-KA2D_features_phys_mean[ind]).div(var)
		ind+=1
		
	# remove features with zero variance and length scale > 5.75 
	length=KA2D_features_phys.shape[1]
	for ind in range(length-1,-1,-1):
		var=KA2D_features_phys_var[ind]
		#print(var)
		if var < 0.00001 :
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

	# create folder for learned models
	PATH='./models_mlp'
	isExist = os.path.exists(PATH)
	if not isExist:
	  os.makedirs(PATH)
	  
	# create folder for results
	PATHR='./results_mlp'
	isExist = os.path.exists(PATHR)
	if not isExist:
	  os.makedirs(PATHR)

	# learn different dynamical descriptors and time scales
	for d in dyn_array:
	
		for t in range(Tstart, Tstop):
						
			print ("Analyze type ", Type, " for observable ", d, " at time ", t)
		
			# choose labels
			KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]
			if d == "BB" :
				KA2D_labels_mod = 1.0 - KA2D_labels
			if d == "MD" :
				KA2D_labels_mod = (tf.tanh(  (KA2D_labels-0.4)*20.0  ) + 1.0)/2.0

			# extract features
			KA2D_features = KA2D_features_phys
			KA2D_features, KA2D_features_test = np.split(KA2D_features, [NSTrain*NLOC])
			
			# calc distance matrices, norms and true label results	
			for s in range(NSTrain + NSTest):
				
				# calc distances
				positions = np.stack((XPOS[NLOC*s:NLOC*(s+1)],YPOS[NLOC*s:NLOC*(s+1)]),axis=1 )
				dist_nd_sq = np.zeros(NLOC * (NLOC - 1) // 2)  # to match the result of pdist
				for dim in range(2):
					pos_1d = positions[:, dim][:, np.newaxis]  # shape (N, 1)
					dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
					dist_1d[dist_1d > boxL * 0.5] -= boxL
					dist_1d[dist_1d < -boxL * 0.5] += boxL
					dist_nd_sq += dist_1d ** 2  # d^2 = dx^2 + dy^2 + dz^2
				dist_nd = squareform(np.sqrt(dist_nd_sq))
				
				# create density matrices
				dist_cg_1 = np.exp(-(dist_nd-2.0)*(dist_nd-2.0)*0.5)
				dist_cg_2 = np.exp(-(dist_nd-5.0)*(dist_nd-5.0)*0.5)
				dist_cg_3 = np.exp(-(dist_nd-8.0)*(dist_nd-8.0)*0.5)
				
				# calc norms
				norm1 = np.sum(dist_cg_1,axis=1)
				norm2 = np.sum(dist_cg_2,axis=1)
				norm3 = np.sum(dist_cg_3,axis=1)
				
				# calc true label result
				KA2D_labels_loc = tf.cast(KA2D_labels_mod[s*NLOC:(s+1)*NLOC], tf.float32)
				# calc true label result
				Mean = np.mean(KA2D_labels_loc)
				Fluct = tf.math.multiply(tf.cast(KA2D_labels_loc-Mean, tf.float32),tf.cast(KA2D_labels_loc-Mean, tf.float32))
				Fluct = tf.reduce_sum(Fluct)
				#KA2D_labels_mean = 0.0
				res1=tf.linalg.matvec(tf.cast(dist_cg_1, tf.float32),KA2D_labels_loc-Mean)
				res2=tf.linalg.matvec(tf.cast(dist_cg_2, tf.float32),KA2D_labels_loc-Mean)
				res3=tf.linalg.matvec(tf.cast(dist_cg_3, tf.float32),KA2D_labels_loc-Mean)
				res1 = tf.math.multiply(res1,KA2D_labels_loc-Mean)
				res2 = tf.math.multiply(res2,KA2D_labels_loc-Mean)
				res3 = tf.math.multiply(res3,KA2D_labels_loc-Mean)
				res1 /= norm1
				res2 /= norm2
				res3 /= norm3
				res1 = tf.reduce_sum(res1)
				res2 = tf.reduce_sum(res2)
				res3 = tf.reduce_sum(res3)
				res1/=Fluct
				res2/=Fluct
				res3/=Fluct
				
				
				res1 = np.full(NLOC,res1)
				res2 = np.full(NLOC,res2)
				res3 = np.full(NLOC,res3)
				
				if s==0 :
					wlist1=dist_cg_1
					wlist2=dist_cg_2
					wlist3=dist_cg_3
					normlist1=norm1
					normlist2=norm2
					normlist3=norm3
					reslist1=res1
					reslist2=res2
					reslist3=res3
				else :
					wlist1=np.concatenate([wlist1,dist_cg_1])
					wlist2=np.concatenate([wlist2,dist_cg_2])
					wlist3=np.concatenate([wlist3,dist_cg_3])
					normlist1=np.concatenate([normlist1,norm1])
					normlist2=np.concatenate([normlist2,norm2])
					normlist3=np.concatenate([normlist3,norm3])
					reslist1=np.concatenate([reslist1,res1])
					reslist2=np.concatenate([reslist2,res2])
					reslist3=np.concatenate([reslist3,res3])
					
			# include data to labels
			KA2D_labels = np.column_stack([KA2D_labels,reslist1,reslist2,reslist3,normlist1,normlist2,normlist3,wlist1,wlist2,wlist3])
			KA2D_labels, KA2D_labels_test = np.split(KA2D_labels, [NSTrain*NLOC])
			KA2D_labels_mod, KA2D_labels_test_mod = np.split(KA2D_labels_mod, [NSTrain*NLOC])
			print(KA2D_labels.shape)
			
			# learn different models
			Nmod = 8
			for mod in range(Nmod):
			
				print ("Analyze type ", Type, " for observable ", d, " at time ", t, " for ML model ", mod)
			
				# define ML model
				Nbottle = 2
				Ninner = 10
				KA2D_model = tf.keras.Sequential()
				
				length=KA2D_features.shape[1]
				struct_input = tf.keras.Input(shape=(length,),name="struct_input")

				bottleneck_layer = Dense(Nbottle, activation="elu",name="bottleneck")(struct_input)
				
				x = Dense(Ninner, activation="elu",name="inner1")(bottleneck_layer)
				x = Dense(Ninner, activation="elu",name="inner2")(x)

				outputs = Dense(1, activation="linear")(x)
				
				KA2D_model = tf.keras.Model(struct_input,outputs, name="mlp_timescale")
				
				# create custom loss function
				Var_rescale = tf.cast(tf.sqrt(tf.math.reduce_variance(KA2D_labels[:,0], axis=-1)), tf.float32)
				def custom_loss(y_true, y_pred):
				    	# loss1, squared distances between labels
				    	loss1 = tf.reduce_mean(tf.abs(y_true[:,0] - y_pred[:,0]), axis=-1)/Var_rescale
				    	
				    	return loss1 
				
				# compile
				KA2D_model.compile(loss=custom_loss,optimizer = optimizers.Adam(5e-4))
				print(KA2D_model.summary())
				
				
				##################### learn free model ###############################
			
			
				# first learn free model for 100 steps
				history = KA2D_model.fit(
				    KA2D_features,
				    KA2D_labels,
				    batch_size=NLOC,
				    epochs=300,
				    shuffle=False
				)				   
				
				# then add callback to enable early stopping
				callbacks = [EarlyStopping(monitor='val_loss', mode='min', patience=15),ModelCheckpoint(filepath=PATH + "/mlpparams_free_{}{}_{}_{}.dat".format(d,t,Type,mod), monitor='val_loss', save_best_only=True, mode='min')]
				KA2D_model.compile(loss=custom_loss,optimizer = optimizers.Adam(2e-4))
				history = KA2D_model.fit(
				    #[KA2D_features,wlist1,wlist2,wlist3,normlist1,normlist2,normlist3,reslist1,reslist2,reslist3],
				    KA2D_features,
				    KA2D_labels,
				    batch_size=NLOC,
				    epochs=100,
				    shuffle=False,
				    validation_data=(KA2D_features_test,KA2D_labels_test),
				    callbacks=callbacks
				)
				
				KA2D_model = load_model(PATH + "/mlpparams_free_{}{}_{}_{}.dat".format(d,t,Type,mod), custom_objects={'custom_loss': custom_loss })
				   
				# eval model
				y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				corr, _ = pearsonr(y_out, KA2D_labels_test[:,0])
				coef, _ = spearmanr(y_out, KA2D_labels_test[:,0])
				
				print("free ",corr,coef)
				
				# eval and prepare for updated loss function
				# print bottleneck layer
				KA2D_model_bottle = tf.keras.Model(struct_input,bottleneck_layer, name="mlp_bottleneck")
				y_out=KA2D_model_bottle.predict_on_batch(KA2D_features_test)
				plt.scatter(y_out[:,0],y_out[:,1], c=KA2D_labels_test[:,0],s=0.2, alpha=0.5, cmap='coolwarm')
				plt.colorbar()
				plt.title("bottleneck {} t={}".format(d, t))
				plt.xlabel("bottleneck d1")
				plt.ylabel("bottleneck d2")
				plt.savefig(PATHR + "/bottleneck_free_{}{}_{}.png".format(d, t, Type))
				plt.clf()
				
				# calc correlation function
				y_out=KA2D_model.predict_on_batch(KA2D_features)[:,0]
				y_out_test=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				if d == "BB" :
					y_out_mod = 1.0 - y_out
					y_out_test_mod = 1.0 - y_out_test
				if d == "MD" :
					y_out_mod = (tf.tanh(  (y_out-0.4)*20.0  ) + 1.0)/2.0
					y_out_test_mod = (tf.tanh(  (y_out_test-0.4)*20.0  ) + 1.0)/2.0
				corr1 = 0.0
				corr2 = 0.0
				corr3 = 0.0
				corr1_true = 0.0
				corr2_true = 0.0
				corr3_true = 0.0
				corr1_test = 0.0
				corr2_test = 0.0
				corr3_test = 0.0
				corr1_true_test = 0.0
				corr2_true_test = 0.0
				corr3_true_test = 0.0
				norm = 0.0
				for s in range(NSTrain + NSTest):
					# calculate cg labels
					if s < NSTrain :
						y_out_loc = tf.cast(y_out_mod[s*NLOC:(s+1)*NLOC], tf.float32)
						KA2D_labels_loc = tf.cast(KA2D_labels_mod[s*NLOC:(s+1)*NLOC], tf.float32)
					else :
						y_out_loc = tf.cast(y_out_test_mod[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC], tf.float32)
						KA2D_labels_loc = tf.cast(KA2D_labels_test_mod[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC], tf.float32)
						#print(KA2D_labels_loc)
					Mean = np.mean(y_out_loc)
					Fluct = tf.math.multiply(tf.cast(y_out_loc-Mean, tf.float32),tf.cast(y_out_loc-Mean, tf.float32))
					Fluct = tf.reduce_sum(Fluct)
					res1=tf.linalg.matvec(tf.cast(wlist1[s*NLOC:(s+1)*NLOC], tf.float32),y_out_loc-Mean)/normlist1[s*NLOC:(s+1)*NLOC]
					res2=tf.linalg.matvec(tf.cast(wlist2[s*NLOC:(s+1)*NLOC], tf.float32),y_out_loc-Mean)/normlist2[s*NLOC:(s+1)*NLOC]
					res3=tf.linalg.matvec(tf.cast(wlist3[s*NLOC:(s+1)*NLOC], tf.float32),y_out_loc-Mean)/normlist3[s*NLOC:(s+1)*NLOC]

					
					# calc correlation function
					res1 = tf.math.multiply(res1,y_out_loc-Mean)
					res2 = tf.math.multiply(res2,y_out_loc-Mean)
					res3 = tf.math.multiply(res3,y_out_loc-Mean)
					res1 = tf.reduce_sum(res1)
					res2 = tf.reduce_sum(res2)
					res3 = tf.reduce_sum(res3)
					
					#res1/=KA2D_labels_mean*KA2D_labels_mean
					#res2/=KA2D_labels_mean*KA2D_labels_mean
					#res3/=KA2D_labels_mean*KA2D_labels_mean
					
					res1/=Fluct
					res2/=Fluct
					res3/=Fluct
					
					#print (res1)
					
					if s < NSTrain :
						corr1 += res1
						corr2 += res2
						corr3 += res3
						corr1_true += reslist1[s*NLOC]
						corr2_true += reslist2[s*NLOC]
						corr3_true += reslist3[s*NLOC]
					else :
						corr1_test += res1
						corr2_test += res2
						corr3_test += res3
						corr1_true_test += reslist1[s*NLOC]
						corr2_true_test += reslist2[s*NLOC]
						corr3_true_test += reslist3[s*NLOC]
				
				print("Correlation Training Set")	
				print (corr1/150, corr1_true/150)
				print (corr2/150, corr2_true/150)
				print (corr3/150, corr3_true/150)
				print("Correlation Test Set")	
				print (corr1_test/150, corr1_true_test/150)
				print (corr2_test/150, corr2_true_test/150)
				print (corr3_test/150, corr3_true_test/150)
					
				print("Histogram")	
				for l in np.linspace(0.1, 1.0, num=9):
					l_test = l;
					
					hist_loc_true = tf.math.square(KA2D_labels[:,0] - l_test)
					hist_loc_true = -hist_loc_true / 0.02
					hist_loc_true = tf.math.exp(hist_loc_true)
					hist_loc_true = tf.reduce_sum(hist_loc_true)
					
					hist_loc_pred = tf.math.square(y_out - l_test)
					hist_loc_pred = -hist_loc_pred / 0.02
					hist_loc_pred = tf.math.exp(hist_loc_pred)
					hist_loc_pred = tf.reduce_sum(hist_loc_pred)
					
					print("l_test:", l_test, hist_loc_pred, hist_loc_true)
					
				print("Var")
				print(np.var(y_out),np.var(KA2D_labels[:,0]))	
					
				print("Mean")
				print(np.mean(y_out),np.mean(KA2D_labels[:,0]))
				
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
				#print((var-mean*mean)/NLOC,var_single-mean_single*mean_single)
				chi4_true=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
				chi4_true_standard=(var-mean*mean)/NLOC
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
				chi4_pred=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
				chi4_pred_standard=(var-mean*mean)/NLOC
				
				print ("chi4, ", chi4_true, chi4_pred, chi4_true_standard, chi4_pred_standard)
				
				
				################################### learn correlation model #####################################

				
				# define updated loss function
				def loss_carrier(factorP, factorC,factorV):
					def loss(y_true, y_pred):
					
						Mean_true = tf.reduce_mean(y_true[:,0], axis=-1)
						Var_true = tf.math.reduce_variance(y_true[:,0], axis=-1)
					
					    	# loss1, squared distances between labels
						loss1 = tf.reduce_mean(tf.abs(y_true[:,0] - y_pred[:,0]), axis=-1)/Var_rescale
						#tf.print("Loss1:", loss1, output_stream=sys.stdout)
						
						# calculate cg labels
						if d == "BB" :
							y_pred_mod = 1.0 - y_pred[:,0]
							y_true_mod = 1.0 - y_true[:,0]
						if d == "MD" :
							y_pred_mod = (tf.tanh(  (y_pred[:,0]-0.4)*20.0  ) + 1.0)/2.0
							y_true_mod = (tf.tanh(  (y_true[:,0]-0.4)*20.0  ) + 1.0)/2.0
						
						Mean = tf.reduce_mean(y_pred_mod, axis=-1)
						Fluct = tf.math.multiply(tf.cast(y_pred_mod-Mean, tf.float32),tf.cast(y_pred_mod-Mean, tf.float32))
						Fluct = tf.reduce_sum(Fluct)
						res1=tf.linalg.matvec(y_true[:,7:NLOC+7],y_pred_mod-Mean)
						res1/=y_true[:,4]
						res2=tf.linalg.matvec(y_true[:,7+NLOC:2*NLOC+7],y_pred_mod-Mean)
						res2/=y_true[:,5]
						res3=tf.linalg.matvec(y_true[:,7+2*NLOC:3*NLOC+7],y_pred_mod-Mean)
						res3/=y_true[:,6]
						
						# calc correlation function
						res1 = tf.math.multiply(res1,y_pred_mod-Mean)
						res2 = tf.math.multiply(res2,y_pred_mod-Mean)
						res3 = tf.math.multiply(res3,y_pred_mod-Mean)
						res1 = tf.reduce_sum(res1)
						res2 = tf.reduce_sum(res2)
						res3 = tf.reduce_sum(res3)
						
						#res1/=Mean_true_mod*Mean_true_mod
						#res2/=Mean_true_mod*Mean_true_mod
						#res3/=Mean_true_mod*Mean_true_mod
						
						res1/=Fluct
						res2/=Fluct
						res3/=Fluct
						
						# loss2, squared distance between cg labels
						loss2_1 = tf.abs((y_true[0,1] - res1)/(y_true[0,1]-y_true[0,3]))
						loss2_2 = tf.abs((y_true[0,2] - res2)/(y_true[0,1]-y_true[0,3]))
						loss2_3 = tf.abs((y_true[0,3]-res3)/(y_true[0,1]-y_true[0,3]))
						#tf.print("Loss2:", loss2_1, loss2_2, loss2_3, output_stream=sys.stdout)
						
						# loss3, histograms
						#loss3 = 0.0
						#for l in np.linspace(0.1, 1.0, num=9):
						#	l_test = l;
							
						#	hist_loc_true = tf.math.square(y_true[:,0] - l_test)
						#	hist_loc_true = -hist_loc_true / 0.02
						#	hist_loc_true = tf.math.exp(hist_loc_true)
						#	hist_loc_true = tf.reduce_sum(hist_loc_true)
							
						#	hist_loc_pred = tf.math.square(y_pred[:,0] - l_test)
						#	hist_loc_pred = -hist_loc_pred / 0.02
						#	hist_loc_pred = tf.math.exp(hist_loc_pred)
						#	hist_loc_pred = tf.reduce_sum(hist_loc_pred)
							
						#	if hist_loc_true > 0.01 :
						#		loss3 += tf.abs(hist_loc_pred - hist_loc_true)/hist_loc_true
								#tf.print("l_test:", l_test,hist_loc_pred,hist_loc_true, tf.abs(hist_loc_pred - hist_loc_true)/hist_loc_true, output_stream=sys.stdout)
						#tf.print("Loss3:", loss3, output_stream=sys.stdout)	
								
						# loss 4 variance
						Var_pred = tf.math.reduce_variance(y_pred[:,0], axis=-1)
						Var_true = tf.math.reduce_variance(y_true[:,0], axis=-1)
						#tf.print("Loss:",Var_pred,Var_true, output_stream=sys.stdout)
						loss4 = tf.abs(Var_pred - Var_true)/Var_true
						#tf.print("Loss4:", loss4, output_stream=sys.stdout)
						
						#loss 5 mean
						#Mean_pred = tf.reduce_mean(y_pred[:,0], axis=-1)
						#tf.print("Loss:",Var_pred,Var_true, output_stream=sys.stdout)
						#loss5 = tf.abs(Mean_pred - Mean_true)/(1.0-Mean_true)
						#tf.print("Loss5:", loss5, output_stream=sys.stdout)
						
						loss = loss1*factorP + (loss2_1 + loss2_2 + loss2_3)*factorC + factorV*loss4 
						#loss = tf.reduce_mean(squared_difference, axis=-1)
						
						return loss
					return loss
				
				# fit with updated loss function
				KA2D_model.compile(loss=loss_carrier(1.0,0.5,1.0),optimizer = optimizers.Adam(4e-5))
				
				# first run without callback
				history = KA2D_model.fit(
				    KA2D_features,
				    KA2D_labels,
				    batch_size=NLOC,
				    epochs=50,
				    shuffle=False
				)
				
				# then include callback
				callbacks = [EarlyStopping(monitor='val_loss', mode='min', patience=15),ModelCheckpoint(filepath=PATH + "/mlpparams_{}{}_{}_{}.dat".format(d,t,Type,mod), monitor='val_loss', save_best_only=True, mode='min')]
				KA2D_model.compile(loss=loss_carrier(1.0,0.5,1.0),optimizer = optimizers.Adam(2e-5))
				history = KA2D_model.fit(
				    KA2D_features,
				    KA2D_labels,
				    batch_size=NLOC,
				    epochs=100,
				    shuffle=False,
				    validation_data=(KA2D_features_test,KA2D_labels_test),
				    callbacks=callbacks
				)
				
				KA2D_model = load_model(PATH + "/mlpparams_{}{}_{}_{}.dat".format( d,t,Type, mod), custom_objects={'loss': loss_carrier(1.0,2.0,2.0)})
				   
				# eval model
				y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				corr, _ = pearsonr(y_out, KA2D_labels_test[:,0])
				coef, _ = spearmanr(y_out, KA2D_labels_test[:,0])
				
				print("correlate ",corr,coef)
				
				# eval and prepare for updated loss function
				# print bottleneck layer
				KA2D_model_bottle = tf.keras.Model(struct_input,bottleneck_layer, name="mlp_bottleneck")
				y_out=KA2D_model_bottle.predict_on_batch(KA2D_features_test)
				plt.scatter(y_out[:,0],y_out[:,1], c=KA2D_labels_test[:,0],s=0.2, alpha=0.5, cmap='coolwarm')
				plt.colorbar()
				plt.title("bottleneck {} t={}".format(d, t))
				plt.xlabel("bottleneck d1")
				plt.ylabel("bottleneck d2")
				plt.savefig(PATHR + "/bottleneck_free_{}{}_{}.png".format(d, t, Type))
				plt.clf()
				
				# calc correlation function
				y_out=KA2D_model.predict_on_batch(KA2D_features)[:,0]
				y_out_test=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				if d == "BB" :
					y_out_mod = 1.0 - y_out
					y_out_test_mod = 1.0 - y_out_test
				if d == "MD" :
					y_out_mod = (tf.tanh(  (y_out-0.4)*20.0  ) + 1.0)/2.0
					y_out_test_mod = (tf.tanh(  (y_out_test-0.4)*20.0  ) + 1.0)/2.0
				corr1 = 0.0
				corr2 = 0.0
				corr3 = 0.0
				corr1_true = 0.0
				corr2_true = 0.0
				corr3_true = 0.0
				corr1_test = 0.0
				corr2_test = 0.0
				corr3_test = 0.0
				corr1_true_test = 0.0
				corr2_true_test = 0.0
				corr3_true_test = 0.0
				norm = 0.0
				for s in range(NSTrain + NSTest):
				
					if s < NSTrain :
						y_out_loc = tf.cast(y_out_mod[s*NLOC:(s+1)*NLOC], tf.float32)
						KA2D_labels_loc = tf.cast(KA2D_labels_mod[s*NLOC:(s+1)*NLOC], tf.float32)
					else :
						y_out_loc = tf.cast(y_out_test_mod[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC], tf.float32)
						KA2D_labels_loc = tf.cast(KA2D_labels_test_mod[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC], tf.float32)
						#print(KA2D_labels_loc)
					Mean = np.mean(y_out_loc)
					Fluct = tf.math.multiply(tf.cast(y_out_loc-Mean, tf.float32),tf.cast(y_out_loc-Mean, tf.float32))
					Fluct = tf.reduce_sum(Fluct)
					res1=tf.linalg.matvec(tf.cast(wlist1[s*NLOC:(s+1)*NLOC], tf.float32),y_out_loc-Mean)/normlist1[s*NLOC:(s+1)*NLOC]
					res2=tf.linalg.matvec(tf.cast(wlist2[s*NLOC:(s+1)*NLOC], tf.float32),y_out_loc-Mean)/normlist2[s*NLOC:(s+1)*NLOC]
					res3=tf.linalg.matvec(tf.cast(wlist3[s*NLOC:(s+1)*NLOC], tf.float32),y_out_loc-Mean)/normlist3[s*NLOC:(s+1)*NLOC]

					
					# calc correlation function
					res1 = tf.math.multiply(res1,y_out_loc-Mean)
					res2 = tf.math.multiply(res2,y_out_loc-Mean)
					res3 = tf.math.multiply(res3,y_out_loc-Mean)
					res1 = tf.reduce_sum(res1)
					res2 = tf.reduce_sum(res2)
					res3 = tf.reduce_sum(res3)
					
					#res1/=KA2D_labels_mean*KA2D_labels_mean
					#res2/=KA2D_labels_mean*KA2D_labels_mean
					#res3/=KA2D_labels_mean*KA2D_labels_mean
					
					res1/=Fluct
					res2/=Fluct
					res3/=Fluct
					
					if s < NSTrain :
						corr1 += res1
						corr2 += res2
						corr3 += res3
						corr1_true += reslist1[s*NLOC]
						corr2_true += reslist2[s*NLOC]
						corr3_true += reslist3[s*NLOC]
					else :
						corr1_test += res1
						corr2_test += res2
						corr3_test += res3
						corr1_true_test += reslist1[s*NLOC]
						corr2_true_test += reslist2[s*NLOC]
						corr3_true_test += reslist3[s*NLOC]
				
				print("Correlation Training Set")	
				print (corr1/150, corr1_true/150)
				print (corr2/150, corr2_true/150)
				print (corr3/150, corr3_true/150)
				print("Correlation Test Set")	
				print (corr1_test/150, corr1_true_test/150)
				print (corr2_test/150, corr2_true_test/150)
				print (corr3_test/150, corr3_true_test/150)
					
				print("Histogram")	
				for l in np.linspace(0.1, 1.0, num=9):
					l_test = l;
					
					hist_loc_true = tf.math.square(KA2D_labels[:,0] - l_test)
					hist_loc_true = -hist_loc_true / 0.02
					hist_loc_true = tf.math.exp(hist_loc_true)
					hist_loc_true = tf.reduce_sum(hist_loc_true)
					
					hist_loc_pred = tf.math.square(y_out - l_test)
					hist_loc_pred = -hist_loc_pred / 0.02
					hist_loc_pred = tf.math.exp(hist_loc_pred)
					hist_loc_pred = tf.reduce_sum(hist_loc_pred)
					
					print("l_test:", l_test, hist_loc_pred, hist_loc_true)
					
				print("Var")
				print(np.var(y_out),np.var(KA2D_labels[:,0]))	
					
				print("Mean")
				print(np.mean(y_out),np.mean(KA2D_labels[:,0]))
				
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
				#print((var-mean*mean)/NLOC,var_single-mean_single*mean_single)
				chi4_true=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
				chi4_true_standard=(var-mean*mean)/NLOC
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
				chi4_pred=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
				chi4_pred_standard=(var-mean*mean)/NLOC
				
				print ("chi4, ", chi4_true, chi4_pred, chi4_true_standard, chi4_pred_standard)		

# execute for the three different types					
for type in range(0,3):
	train_model(NArray[type],NameArray[type])

