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

N=4096
NType=2
NArray = [3277, 819]
NameArray = ["type0", "type1"]
boxL=15.041

Temp=0.44
NT=10
NSTrain=30
NSTest=30
Tstart=1
Tstop=NT

tanhSTEP = 0.44
tanhWIDTH = 20.0

dyn_array = ['MD']

def train_model(NLOC,Type):

	# read simulation data
	# labels
	KA_labels_full = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T"+ str(Temp) +"_bapst/ml_labels_" + Type + ".csv",nrows=(NSTrain+NSTest)*NLOC+2)
	KA_labels_full = KA_labels_full.iloc[: , :-1]
	KA_labels_full_mean = KA_labels_full.values[0]
	KA_labels_full_var = KA_labels_full.values[1]
	KA_labels_full.drop([0,1], inplace=True)

	# features
	KA_features_phys = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T" + str(Temp) +"_bapst/ml_struct_" + Type + ".csv",nrows=(NSTrain+NSTest)*NLOC+3)
	KA_features_phys = KA_features_phys.iloc[: , :-1]

	# extract positions
	XPOS = KA_features_phys.pop('NDim0')
	XPOS.drop([0,1,2], inplace=True)
	YPOS = KA_features_phys.pop('NDim1')
	YPOS.drop([0,1,2], inplace=True)
	ZPOS = KA_features_phys.pop('NDim2')
	ZPOS.drop([0,1,2], inplace=True)
	KA_features_phys_length = KA_features_phys.values[0]
	KA_features_phys_mean = KA_features_phys.values[1]
	KA_features_phys_var = KA_features_phys.values[2]
	KA_features_phys.drop([0,1,2], inplace=True)
	
	# rescale features
	labels=KA_features_phys.columns.values
	ind=0
	for var in KA_features_phys_var:
		#print(ind)
		KA_features_phys[labels[ind]] = KA_features_phys[labels[ind]].add(-KA_features_phys_mean[ind]).div(var)
		ind+=1
		
	# remove features with zero variance and length scale > 5.75 
	length=KA_features_phys.shape[1]
	for ind in range(length-1,-1,-1):
		var=KA_features_phys_var[ind]
		#print(var)
		if var < 0.00001 :
			KA_features_phys.drop(KA_features_phys.columns[ind], axis=1, inplace=True)
			KA_features_phys_length=np.delete(KA_features_phys_length,ind)
		ind+=1

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
			KA_labels = KA_labels_full["{}{}".format(d, t)]
			#KA_labels_Trans = np.where(KA_labels < 0.44, 0.0, 1.0)
			KA_labels_tanh = (tf.tanh(  (KA_labels-tanhSTEP)*tanhWIDTH  ) + 1.0)/2.0
			
			#print(KA_labels_Trans)
			print(KA_labels_tanh)

			# extract features
			KA_features = KA_features_phys
			KA_features, KA_features_test = np.split(KA_features, [NSTrain*NLOC])
			KA_features.index -= 3
			KA_features_test.index -= 3
			
			# calc distance matrices, norms and true label results	
			for s in range(NSTrain + NSTest):
				print(s)
				# calc distances
				positions = np.stack((XPOS[NLOC*s:NLOC*(s+1)],YPOS[NLOC*s:NLOC*(s+1)],ZPOS[NLOC*s:NLOC*(s+1)]),axis=1 )
				dist_nd_sq = np.zeros(NLOC * (NLOC - 1) // 2)  # to match the result of pdist
				for dim in range(3):
					pos_1d = positions[:, dim][:, np.newaxis]  # shape (N, 1)
					dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
					dist_1d[dist_1d > boxL * 0.5] -= boxL
					dist_1d[dist_1d < -boxL * 0.5] += boxL
					dist_nd_sq += dist_1d ** 2  # d^2 = dx^2 + dy^2 + dz^2
				dist_nd = squareform(np.sqrt(dist_nd_sq))
				
				# create density matrices
				dist_cg_1 = np.exp(-(dist_nd-2.0)*(dist_nd-2.0)*0.5)
				dist_cg_2 = np.exp(-(dist_nd-4.0)*(dist_nd-4.0)*0.5)
				dist_cg_3 = np.exp(-(dist_nd-6.0)*(dist_nd-6.0)*0.5)
				
				# calc norms
				norm1 = np.sum(dist_cg_1,axis=1)
				norm2 = np.sum(dist_cg_2,axis=1)
				norm3 = np.sum(dist_cg_3,axis=1)
				
				# calc true label result
				KA_labels_tanh_loc = tf.cast(KA_labels_tanh[s*NLOC:(s+1)*NLOC], tf.float32)
				Mean = np.mean(KA_labels_tanh_loc)
				Fluct = tf.math.multiply(tf.cast(KA_labels_tanh_loc-Mean, tf.float32),tf.cast(KA_labels_tanh_loc-Mean, tf.float32))
				Fluct = tf.reduce_sum(Fluct)
				res1=tf.linalg.matvec(tf.cast(dist_cg_1, tf.float32),KA_labels_tanh_loc-Mean)
				res2=tf.linalg.matvec(tf.cast(dist_cg_2, tf.float32),KA_labels_tanh_loc-Mean)
				res3=tf.linalg.matvec(tf.cast(dist_cg_3, tf.float32),KA_labels_tanh_loc-Mean)
				res1 = tf.math.multiply(res1,KA_labels_tanh_loc-Mean)
				res2 = tf.math.multiply(res2,KA_labels_tanh_loc-Mean)
				res3 = tf.math.multiply(res3,KA_labels_tanh_loc-Mean)
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
					
			KA_labels, KA_labels_test = np.split(KA_labels, [NSTrain*NLOC])
			KA_labels_tanh, KA_labels_tanh_test = np.split(KA_labels_tanh, [NSTrain*NLOC])
			
			# create data generator
			"""class CustomDataGen(tf.keras.utils.Sequence):
			
				def __init__(self, df_input, df_label,batch_size):
					self.df_input = df_input.copy()
					self.df_label = df_label.copy()
					self.batch_size = batch_size
					self.n = len(self.df_input)
					
					
				def __getitem__(self, index):
					#tf.print("Test Generate", output_stream=sys.stdout)
					batch_input = self.df_input[index * self.batch_size:(index + 1) * self.batch_size]
					#tf.print("Test Generate", batch_input.shape,batch_mean_input.shape, batch_label.shape, output_stream=sys.stdout)
					batch_label = self.df_label[index * self.batch_size:(index + 1) * self.batch_size]
					#ind_loc = batch_label.pop(7)
					#tf.print("Test Generate",ind_loc)
					#wlist1_loc = wlist1[int(ind_loc[0]):int(ind_loc[self.batch_size-1])+1]
					#wlist2_loc = wlist2[int(ind_loc[0]):int(ind_loc[self.batch_size-1])+1]
					#wlist3_loc = wlist3[int(ind_loc[0]):int(ind_loc[self.batch_size-1])+1]
					#batch_label = pd.concat([batch_label, pd.DataFrame(wlist1_loc), pd.DataFrame(wlist2_loc), pd.DataFrame(wlist3_loc)], axis=1)
					#tf.print("Test Generate",batch_label.shape)
					#tf.print("Test Generate", batch_input.shape,batch_mean_input.shape,batch_temp_input.shape, batch_label.shape, output_stream=sys.stdout)
					return batch_input, batch_label
					
				def __len__(self):
					return self.n // self.batch_size"""
					
			#traingen = CustomDataGen(KA_features,KA_labels,batch_size=NLOC)
			#valgen = CustomDataGen(KA_features_test,KA_labels_test,batch_size=NLOC)
			
			# define ML model
			Nbottle = 2
			Ninner = 10
			KA_model = tf.keras.Sequential()
			
			length=KA_features.shape[1]
			struct_input = tf.keras.Input(shape=(length,),name="struct_input")

			bottleneck_layer = Dense(Nbottle, activation="elu",name="bottleneck")(struct_input)
			
			x = Dense(Ninner, activation="elu",name="inner1")(bottleneck_layer)
			x = Dense(Ninner, activation="elu",name="inner2")(x)

			outputs = Dense(1, activation="linear")(x)
			
			KA_model = tf.keras.Model(struct_input,outputs, name="mlp_timescale")
			
			# create custom loss function
			Var_rescale = tf.cast(tf.sqrt(tf.math.reduce_variance(KA_labels.to_numpy(), axis=-1)), tf.float32)
			def custom_loss(y_true, y_pred):
			    	loss1 = tf.reduce_mean(tf.abs(y_true - y_pred), axis=-1)/Var_rescale
			    	return loss1 
			    	
			
			
			# compile
			KA_model.compile(loss=custom_loss,optimizer = optimizers.Adam(5e-4))
			print(KA_model.summary())
			
			
			##################### learn free model ###############################
			
			
			# first learn free model for 100 steps
			history = KA_model.fit(
			    KA_features,
			    KA_labels,
			    batch_size=NLOC,
			    epochs=300,
			    shuffle=False
			)				   
			
			# then add callback to enable early stopping
			callbacks = [EarlyStopping(monitor='val_loss', mode='min', patience=15),ModelCheckpoint(filepath=PATH + "/mlpparams_free_{}{}_{}.dat".format(d,t,Type), monitor='val_loss', save_best_only=True, mode='min')]
			KA_model.compile(loss=custom_loss,optimizer = optimizers.Adam(2e-4))
			history = KA_model.fit(
			    KA_features,
			    KA_labels,
			    batch_size=NLOC,
			    epochs=100,
			    validation_data=(KA_features_test,KA_labels_test),
			    callbacks=callbacks
			)
			
			KA_model = load_model(PATH + "/mlpparams_free_{}{}_{}.dat".format(d,t,Type), custom_objects={'custom_loss': custom_loss })
			   
			# eval model
			y_out=KA_model.predict_on_batch(KA_features_test)[:,0]
			y_out_tanh_test = (tf.tanh(  (y_out-tanhSTEP)*tanhWIDTH  ) + 1.0)/2.0
			corr, _ = pearsonr(y_out, KA_labels_test)
			coef, _ = spearmanr(y_out, KA_labels_test)
			
			print("free ",corr,coef)
			
			# eval and prepare for updated loss function
			# print bottleneck layer
			KA_model_bottle = tf.keras.Model(struct_input,bottleneck_layer, name="mlp_bottleneck")
			y_out=KA_model_bottle.predict_on_batch(KA_features_test)
			plt.scatter(y_out[:,0],y_out[:,1], c=KA_labels_test,s=0.2, alpha=0.5, cmap='coolwarm')
			plt.colorbar()
			plt.title("bottleneck {} t={}".format(d, t))
			plt.xlabel("bottleneck d1")
			plt.ylabel("bottleneck d2")
			plt.savefig(PATHR + "/bottleneck_free_{}{}_{}.png".format(d, t, Type))
			plt.clf()
			
			# calc correlation function
			y_out=KA_model.predict_on_batch(KA_features)[:,0]
			y_out_tanh = (tf.tanh(  (y_out-tanhSTEP)*tanhWIDTH  ) + 1.0)/2.0
			y_out_test=KA_model.predict_on_batch(KA_features_test)[:,0]
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
					y_out_loc = tf.cast(y_out_tanh[s*NLOC:(s+1)*NLOC], tf.float32)
				else :
					y_out_loc = tf.cast(y_out_tanh_test[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC], tf.float32)

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
			
			print("Correlation Training Set PRED TRUE")	
			print (corr1/150, corr1_true/150)
			print (corr2/150, corr2_true/150)
			print (corr3/150, corr3_true/150)
			print("Correlation Test Set PRED TRUE")	
			print (corr1_test/150, corr1_true_test/150)
			print (corr2_test/150, corr2_true_test/150)
			print (corr3_test/150, corr3_true_test/150)
				
			print("Histogram")	
			for l in np.linspace(0.1, 1.0, num=9):
				l_test = l;
				
				hist_loc_true = tf.math.square(KA_labels - l_test)
				hist_loc_true = -hist_loc_true / 0.02
				hist_loc_true = tf.math.exp(hist_loc_true)
				hist_loc_true = tf.reduce_sum(hist_loc_true)
				
				hist_loc_pred = tf.math.square(y_out - l_test)
				hist_loc_pred = -hist_loc_pred / 0.02
				hist_loc_pred = tf.math.exp(hist_loc_pred)
				hist_loc_pred = tf.reduce_sum(hist_loc_pred)
				
				print("l_test:", l_test, hist_loc_pred, hist_loc_true)
				
			print("Var")
			print(np.var(y_out),np.var(KA_labels))
			
			print("Var Tanh")
			print(np.var(y_out_tanh_test),np.var(KA_labels_tanh_test))	
				
			print("Mean")
			print(np.mean(y_out),np.mean(KA_labels))
			
			print("Mean Tanh")
			print(np.mean(y_out_tanh_test),np.mean(KA_labels_tanh_test))
			
			# calc chi_4
			var = mean = 0
			var_single = mean_single = 0
			for s in range(NSTest):
				qtot_true = 0
				labels = KA_labels_tanh_test[s*NLOC:(s+1)*NLOC]
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
				labels = y_out_tanh_test[s*NLOC:(s+1)*NLOC]
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
			
								
			# include data to labels
			KA_labels = np.concatenate( [KA_labels,  KA_labels_test])
			KA_labels = np.column_stack([KA_labels,reslist1,reslist2,reslist3,normlist1,normlist2,normlist3])
			KA_labels, KA_labels_test = np.split(KA_labels, [NSTrain*NLOC])
			#print(KA_labels.shape)

			
			# define updated loss function
			def loss_carrier(factorP, factorC,factorV):
				def loss(y_true, y_pred):
				
					Mean_true = tf.reduce_mean(y_true[:,0], axis=-1)
					Var_true = tf.math.reduce_variance(y_true[:,0], axis=-1)
				
				    	# loss1, squared distances between labels
					loss1 = tf.reduce_mean(tf.abs(y_true[:,0] - y_pred[:,0]), axis=-1)/Var_rescale
					#tf.print("Loss1:", loss1, output_stream=sys.stdout)
					
					# calculate cg labels
					y_pred_mod = (tf.tanh(  (y_pred[:,0]-0.44)*20.0  ) + 1.0)/2.0
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
					loss2_3 = tf.abs((y_true[0,3] - res3)/(y_true[0,1]-y_true[0,3]))
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
					#y_true_tanh = (tf.tanh(  (y_true[:,0]-tanhSTEP)*tanhWIDTH  ) + 1.0)/2.0
					#Var_true_tanh = tf.math.reduce_variance(y_true_tanh, axis=-1)
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
				
				
			# create data generator
			class CustomDataGen2(tf.keras.utils.Sequence):
			
				def __init__(self, df_input, df_label,batch_size):
					self.df_input = df_input.copy()
					self.df_label = df_label.copy()
					self.batch_size = batch_size
					self.n = len(self.df_input)
					
				def __getitem__(self, index):
					#tf.print("Test Generate", output_stream=sys.stdout)
					batch_input = self.df_input[index * self.batch_size:(index + 1) * self.batch_size]

					batch_label = self.df_label[index * self.batch_size:(index + 1) * self.batch_size]
					#ind_loc = batch_label.pop(7)
					#tf.print("Test Generate",ind_loc)
					wlist1_loc = wlist1[index * self.batch_size:(index + 1) * self.batch_size]
					wlist2_loc = wlist2[index * self.batch_size:(index + 1) * self.batch_size]
					wlist3_loc = wlist3[index * self.batch_size:(index + 1) * self.batch_size]
					#ind_loc = None
					batch_label = np.column_stack([batch_label, wlist1_loc, wlist2_loc, wlist3_loc])
					#tf.print("Test Generate",batch_label.shape)
					#tf.print("Test Generate", batch_input.shape, batch_label.shape, output_stream=sys.stdout)
					return np.array(batch_input), np.array(batch_label)
					
				def __len__(self):
					return self.n // self.batch_size
				
			traingen = CustomDataGen2(KA_features,KA_labels,batch_size=NLOC)
			valgen = CustomDataGen2(KA_features_test,KA_labels_test,batch_size=NLOC)
			
			# fit with updated loss function
			KA_model.compile(loss=loss_carrier(1.0,0.5,1.0),optimizer = optimizers.Adam(4e-5))
			
			# first run without callback
			history = KA_model.fit(
			    traingen,
			    epochs=50,
			    shuffle=False
			)
			
			# then include callback
			callbacks = [EarlyStopping(monitor='val_loss', mode='min', patience=15),ModelCheckpoint(filepath=PATH + "/mlpparams_{}{}_{}.dat".format(d,t,Type), monitor='val_loss', save_best_only=True, mode='min')]
			KA_model.compile(loss=loss_carrier(1.0,0.5,1.0),optimizer = optimizers.Adam(2e-5))
			history = KA_model.fit(
			    traingen,
			    epochs=100,
			    shuffle=False,
			    validation_data=valgen,
			    callbacks=callbacks
			)
			
			KA_model = load_model(PATH + "/mlpparams_{}{}_{}.dat".format( d,t,Type), custom_objects={'loss': loss_carrier(1.0,2.0,2.0)})
			   
			# eval model
			y_out=KA_model.predict_on_batch(KA_features_test)[:,0]
			y_out_tanh_test = (tf.tanh(  (y_out-tanhSTEP)*tanhWIDTH  ) + 1.0)/2.0

			corr, _ = pearsonr(y_out, KA_labels_test[:,0])
			coef, _ = spearmanr(y_out, KA_labels_test[:,0])
			
			print("corr ",corr,coef, np.mean(y_out), np.max(y_out[:NLOC]), np.max(KA_labels_test[:NLOC,0]))
			
			# print bottleneck layer
			KA_model_bottle = tf.keras.Model(struct_input,bottleneck_layer, name="mlp_bottleneck")
			y_out=KA_model_bottle.predict_on_batch(KA_features_test)
			plt.scatter(y_out[:,0],y_out[:,1], c=KA_labels_test[:,0],s=0.2, alpha=0.5, cmap='coolwarm')
			plt.colorbar()
			plt.title("bottleneck {} t={}".format(d, t))
			plt.xlabel("bottleneck d1")
			plt.ylabel("bottleneck d2")
			plt.savefig(PATHR + "/bottleneck_{}{}_{}.png".format(d, t, Type))
			plt.clf()
			
			# calc correlation function
			y_out=KA_model.predict_on_batch(KA_features)[:,0]
			print(y_out)
			print (np.max(y_out))
			y_out_tanh = (tf.tanh(  (y_out-tanhSTEP)*tanhWIDTH  ) + 1.0)/2.0
			y_out_test=KA_model.predict_on_batch(KA_features_test)[:,0]
			print(y_out_test)
			print (np.max(y_out_test))
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
					y_out_loc = tf.cast(y_out_tanh[s*NLOC:(s+1)*NLOC], tf.float32)
				else :
					y_out_loc = tf.cast(y_out_tanh_test[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC], tf.float32)

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
			
			print("Correlation Training Set PRED TRUE")	
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
				
				hist_loc_true = tf.math.square(KA_labels[:,0] - l_test)
				hist_loc_true = -hist_loc_true / 0.02
				hist_loc_true = tf.math.exp(hist_loc_true)
				hist_loc_true = tf.reduce_sum(hist_loc_true)
				
				hist_loc_pred = tf.math.square(y_out - l_test)
				hist_loc_pred = -hist_loc_pred / 0.02
				hist_loc_pred = tf.math.exp(hist_loc_pred)
				hist_loc_pred = tf.reduce_sum(hist_loc_pred)
				
				print("l_test:", l_test, hist_loc_pred, hist_loc_true)
				
			print("Var")
			print(np.var(y_out),np.var(KA_labels[:,0]))
			
			print("Var Tanh")
			print(np.var(y_out_tanh_test),np.var(KA_labels_tanh_test))	
				
			print("Mean")
			print(np.mean(y_out),np.mean(KA_labels[:,0]))
			
			print("Mean Tanh")
			print(np.mean(y_out_tanh_test),np.mean(KA_labels_tanh_test))
			
			# calc chi_4
			var = mean = 0
			var_single = mean_single = 0
			for s in range(NSTest):
				qtot_true = 0
				labels = KA_labels_tanh_test[s*NLOC:(s+1)*NLOC]
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
				labels = y_out_tanh_test[s*NLOC:(s+1)*NLOC]
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
			
			KA_model = None
			KA_labels = None
			wlist1 = None
			wlist2 = None
			wlist3 = None

# execute for the three different types					
for type in range(1):
	train_model(NArray[type],NameArray[type])
		

