from signal import SIGUSR1
import pandas as pd
import numpy as np
from numpy import sqrt
from numpy import exp
import os
import sys

# fitting
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge

# Make numpy values easier to read.
np.set_printoptions(precision=3, suppress=True)

# plotting:
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.dpi']=300 # highres display

import tensorflow_probability as tfp
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

Temp=0.251
NT=42
NSTrain=150
NSTest=150

dyn_array = ['BB']
mode_array = ['basic_phys_inh']
#dyn_array = ['MD', 'LOG(FRES)']
#mode_array = ['basic', 'basic_phys_inh', 'basic_phys']


def train_model(NLOC,Type):
	# read simulation data
	KA2D_labels_full = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T"+ str(Temp) +"/ml_labels_" + Type + ".csv",nrows=(NSTrain+NSTest)*NLOC+2)
	KA2D_labels_full = KA2D_labels_full.iloc[: , :-1]
	KA2D_labels_full_mean = KA2D_labels_full.values[0]
	KA2D_labels_full_var = KA2D_labels_full.values[1]
	#print(KA2D_features_phys)
	KA2D_labels_full.drop([0,1], inplace=True)
	#labels=KA2D_labels_full.columns.values
	#ind=0
	#for var in KA2D_labels_full_var:
		#print(ind)
	#	KA2D_labels_full[labels[ind]] = KA2D_labels_full[labels[ind]].add(-KA2D_labels_full_mean[ind]).div(var)
#		ind+=1
#	print(KA2D_labels_full)

	KA2D_features_phys = pd.read_csv("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"/ml_struct_" + Type + ".csv",nrows=(NSTrain+NSTest)*NLOC+3)
	KA2D_features_phys = KA2D_features_phys.iloc[: , :-1]
	KA2D_features_phys=KA2D_features_phys.loc[:,~KA2D_features_phys.columns.str.contains('EPOT_VAR2CGP_ALL')]
	KA2D_features_phys=KA2D_features_phys.loc[:,~KA2D_features_phys.columns.str.contains('PERI_CGP_ALL')]
	XPOS = KA2D_features_phys.pop('DIM0')
	XPOS.drop([0,1,2], inplace=True)
	print (XPOS)
	YPOS = KA2D_features_phys.pop('DIM1')
	YPOS.drop([0,1,2], inplace=True)
	KA2D_features_phys_length = KA2D_features_phys.values[0]
	KA2D_features_phys_mean = KA2D_features_phys.values[1]
	KA2D_features_phys_var = KA2D_features_phys.values[2]
	#print(KA2D_features_phys)
	KA2D_features_phys.drop([0,1,2], inplace=True)
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
	  
	  
	
	# remove all descriptors with length scale > cg
	for cg in np.arange(4.0, 3.5, -1.0):
		length=KA2D_features_phys.shape[1]
		for ind in range(length-1,-1,-1):
			lengthscale=KA2D_features_phys_length[ind]
			if lengthscale > cg :
				KA2D_features_phys.drop(KA2D_features_phys.columns[ind], axis=1, inplace=True)
				KA2D_features_phys_length=np.delete(KA2D_features_phys_length,ind)
			ind+=1
		print(KA2D_features_phys)

		for m in mode_array:
			for d in dyn_array:
				for t in range(25, NT):
					
					print ("Analyze type ", Type, " mode ", m, " for observable ", d, " cg length ", cg," at time ", t)
				
					# choose labels
					KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]
					print(KA2D_labels)

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
						dist_cg_2 = np.exp(-(dist_nd-4.0)*(dist_nd-4.0)*0.5)
						dist_cg_3 = np.exp(-(dist_nd-6.0)*(dist_nd-6.0)*0.5)
						
						# calc norms
						norm1 = np.sum(dist_cg_1,axis=1)
						norm2 = np.sum(dist_cg_2,axis=1)
						norm3 = np.sum(dist_cg_3,axis=1)
						
						# calc true label result
						Mean = np.mean(KA2D_labels[s*NLOC:(s+1)*NLOC])
						Fluct = tf.math.multiply(tf.cast(KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean, tf.float32),tf.cast(KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean, tf.float32))
						Fluct = tf.reduce_sum(Fluct)
						res1=tf.linalg.matvec(tf.cast(dist_cg_1, tf.float32),tf.cast(KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean, tf.float32))
						res2=tf.linalg.matvec(tf.cast(dist_cg_2, tf.float32),tf.cast(KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean, tf.float32))
						res3=tf.linalg.matvec(tf.cast(dist_cg_3, tf.float32),tf.cast(KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean, tf.float32))
						res1 = tf.math.multiply(res1,KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean)
						res2 = tf.math.multiply(res2,KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean)
						res3 = tf.math.multiply(res3,KA2D_labels[s*NLOC:(s+1)*NLOC]-Mean)
						res1 /= norm1
						res2 /= norm2
						res3 /= norm3
						res1 = tf.reduce_sum(res1)
						res2 = tf.reduce_sum(res2)
						res3 = tf.reduce_sum(res3)
						res1 = np.full(NLOC,res1/Fluct)
						res2 = np.full(NLOC,res2/Fluct)
						res3 = np.full(NLOC,res3/Fluct)
						
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
					print(KA2D_labels.shape)
					
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
					def custom_loss(y_true, y_pred):
					    	# loss1, squared distances between labels
						squared_difference = tf.abs(y_true[:,0] - y_pred[:,0])
						loss = tf.reduce_mean(squared_difference, axis=-1) 
						return loss
					
					# compile
					KA2D_model.compile(loss=custom_loss,optimizer = optimizers.Adam(1e-3))
					print(KA2D_model.summary())
					
					
					##################### learn free model ###############################
					
					
					# first learn free model for 100 steps
					history = KA2D_model.fit(
					    KA2D_features,
					    KA2D_labels,
					    batch_size=NLOC,
					    epochs=100,
					    shuffle=False
					)				   
					loss_loc=history.history['loss'][99]
					
					# then add callback to enable early stopping
					# create custom metric function
					def pearson_loss(y_true, y_pred):
						loss=tfp.stats.covariance((y_true[:,0])[:,None],(y_pred[:,0])[:,None])
						norm1= tf.math.sqrt(tfp.stats.covariance((y_true[:,0])[:,None]))
						norm2= tf.math.sqrt(tfp.stats.covariance((y_pred[:,0])[:,None]))
						return loss/(norm1*norm2)
					KA2D_model.compile(loss=custom_loss,optimizer = optimizers.Adam(1e-3), metrics=['accuracy',pearson_loss])
					callbacks = [EarlyStopping(monitor='val_pearson_loss', mode='max', patience=20),ModelCheckpoint(filepath=PATH + "/mlpparams_free_{}{}{}_{}_{}.dat".format(m, d,t,Type,cg), monitor='val_pearson_loss', save_best_only=True, mode='max')]
					history = KA2D_model.fit(
					    #[KA2D_features,wlist1,wlist2,wlist3,normlist1,normlist2,normlist3,reslist1,reslist2,reslist3],
					    KA2D_features,
					    KA2D_labels,
					    batch_size=NLOC,
					    epochs=1000,
					    shuffle=False,
					    validation_data=[KA2D_features_test,KA2D_labels_test],
					    callbacks=callbacks
					)
					
					KA2D_model = load_model(PATH + "/mlpparams_free_{}{}{}_{}_{}.dat".format(m, d,t,Type,cg), custom_objects={'pearson_loss': pearson_loss, 'custom_loss': custom_loss })
					   
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
					plt.title("bottleneck {} {} t={}".format(m, d, t))
					plt.xlabel("bottleneck d1")
					plt.ylabel("bottleneck d2")
					plt.savefig(PATHR + "/bottleneck_free_{}{}{}_{}_{}.png".format(m, d, t, Type,cg))
					plt.clf()
					
					# calc correlation function
					y_out=KA2D_model.predict_on_batch(KA2D_features)
					y_out_test=KA2D_model.predict_on_batch(KA2D_features_test)
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
							y_out_loc = tf.cast(y_out[s*NLOC:(s+1)*NLOC,0], tf.float32)
						else :
							y_out_loc = tf.cast(y_out_test[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC,0], tf.float32)
						Mean = tf.reduce_mean(y_out_loc, axis=-1)
						Fluct = tf.math.multiply(y_out_loc-Mean,y_out_loc-Mean)
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
						res1 = res1/Fluct
						res2 = res2/Fluct
						res3 = res3/Fluct
						
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
						
					print (corr1/150, corr1_true/150)
					print (corr2/150, corr2_true/150)
					print (corr3/150, corr3_true/150)
					print (corr1_test/150, corr1_true_test/150)
					print (corr2_test/150, corr2_true_test/150)
					print (corr3_test/150, corr3_true_test/150)
					diff1 = corr1/150 - corr1_true/150
					diff2 = corr2/150 - corr2_true/150
					diff3 = corr3/150 - corr3_true/150
					if t< 35:
						fac1 = 0.15*loss_loc/diff1
						fac2 = 0.15*loss_loc/diff2
						fac3 = 0.15*loss_loc/diff3
					else :
						fac1 = 0.1*loss_loc/diff1
						fac2 = 0.1*loss_loc/diff2
						fac3 = 0.1*loss_loc/diff3
						
					print ("Change Loss ", fac1,fac2,fac3)
					
					
					################################### learn correlation model #####################################

					
					# define updated loss function
					def loss_carrier(factor1, factor2, factor3, NLOC_loc):
						def loss(y_true, y_pred):
						    	# loss1, squared distances between labels
							squared_difference = tf.abs(y_true[:,0] - y_pred[:,0])
							
							# calculate cg labels
							Mean = tf.reduce_mean(y_pred[:,0], axis=-1)
							Fluct = tf.math.multiply(tf.cast(y_pred[:,0]-Mean, tf.float32),tf.cast(y_pred[:,0]-Mean, tf.float32))
							Fluct = tf.reduce_sum(Fluct)
							#tf.print("tensors:", y_true.get_shape().as_list()[1], NLOC, output_stream=sys.stdout)

							#tf.print("NLOC:", NLOC_loc, output_stream=sys.stdout)
							#tf.print("Shape:", vec1.shape, y_pred[:,0].shape, output_stream=sys.stdout)

							res1=tf.linalg.matvec(y_true[:,7:NLOC_loc+7],(y_pred[:,0])-Mean)
							res1/=y_true[:,4]
							res2=tf.linalg.matvec(y_true[:,7+NLOC_loc:2*NLOC_loc+7],y_pred[:,0]-Mean)
							res2/=y_true[:,5]
							res3=tf.linalg.matvec(y_true[:,7+2*NLOC_loc:3*NLOC_loc+7],y_pred[:,0]-Mean)
							res3/=y_true[:,6]
							
							# calc correlation function
							res1 = tf.math.multiply(res1,y_pred[:,0]-Mean)
							res2 = tf.math.multiply(res2,y_pred[:,0]-Mean)
							res3 = tf.math.multiply(res3,y_pred[:,0]-Mean)
							res1 = tf.reduce_sum(res1)
							res2 = tf.reduce_sum(res2)
							res3 = tf.reduce_sum(res3)
							res1 = res1/Fluct
							res2 = res2/Fluct
							res3 = res3/Fluct
							
							# loss2, squared distance between cg labels
							squared_difference1 = tf.abs(y_true[0,1] - res1)
							squared_difference2 = tf.abs(y_true[0,2] - res2)
							squared_difference3 = tf.abs(y_true[0,3] - res3)
							#tf.print("Loss:", tf.reduce_mean(squared_difference, axis=-1), factor1*squared_difference1,factor2*squared_difference2, factor3*squared_difference3, output_stream=sys.stdout)
							
							loss = tf.reduce_mean(squared_difference, axis=-1) + factor1*squared_difference1 + factor2*squared_difference2 + factor3*squared_difference3
							#loss = tf.reduce_mean(squared_difference, axis=-1)
							
							return loss
						return loss
					
					# fit with updated loss function
					#KA2D_model.compile(loss=loss_carrier(0.0,0.0,0.0, NLOC),optimizer = optimizers.Adam(1e-3), metrics=['accuracy',pearson_loss])
					KA2D_model.compile(loss=loss_carrier(fac1,fac2,fac3, NLOC),optimizer = optimizers.Adam(1e-3), metrics=['accuracy',pearson_loss])
					
					# first run without callback
					history = KA2D_model.fit(
					    KA2D_features,
					    KA2D_labels,
					    batch_size=NLOC,
					    epochs=25,
					    shuffle=False
					)
					
					# then include callback
					callbacks = [EarlyStopping(monitor='val_pearson_loss', mode='max', patience=20),ModelCheckpoint(filepath=PATH + "/mlpparams_{}{}{}_{}_{}.dat".format(m, d,t,Type,cg), monitor='val_pearson_loss', save_best_only=True, mode='max')]
					history = KA2D_model.fit(
					    KA2D_features,
					    KA2D_labels,
					    batch_size=NLOC,
					    epochs=1000,
					    shuffle=False,
					    validation_data=[KA2D_features_test,KA2D_labels_test],
					    callbacks=callbacks
					)
					
					KA2D_model = load_model(PATH + "/mlpparams_{}{}{}_{}_{}.dat".format(m, d,t,Type,cg), custom_objects={'loss': loss_carrier(fac1,fac2,fac3, NLOC), 'pearson_loss': pearson_loss })
					   
					y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
					corr, _ = pearsonr(y_out, KA2D_labels_test[:,0])
					coef, _ = spearmanr(y_out, KA2D_labels_test[:,0])
					
					print("corr ",corr,coef)
					
					# print bottleneck layer
					KA2D_model_bottle = tf.keras.Model(struct_input,bottleneck_layer, name="mlp_bottleneck")
					y_out=KA2D_model_bottle.predict_on_batch(KA2D_features_test)
					plt.scatter(y_out[:,0],y_out[:,1], c=KA2D_labels_test[:,0],s=0.2, alpha=0.5, cmap='coolwarm')
					plt.colorbar()
					plt.title("bottleneck {} {} t={}".format(m, d, t))
					plt.xlabel("bottleneck d1")
					plt.ylabel("bottleneck d2")
					plt.savefig(PATHR + "/bottleneck_{}{}{}_{}_{}.png".format(m, d, t, Type,cg))
					plt.clf()
					
					# calc correlation function
					y_out=KA2D_model.predict_on_batch(KA2D_features)
					y_out_test=KA2D_model.predict_on_batch(KA2D_features_test)
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
							y_out_loc = tf.cast(y_out[s*NLOC:(s+1)*NLOC,0], tf.float32)
						else :
							y_out_loc = tf.cast(y_out_test[(s-NSTrain)*NLOC:(s-NSTrain+1)*NLOC,0], tf.float32)
						Mean = tf.reduce_mean(y_out_loc, axis=-1)
						Fluct = tf.math.multiply(y_out_loc-Mean,y_out_loc-Mean)
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
						res1 = res1/Fluct
						res2 = res2/Fluct
						res3 = res3/Fluct
						
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
						
					print (corr1/150, corr1_true/150)
					print (corr2/150, corr2_true/150)
					print (corr3/150, corr3_true/150)
					print (corr1_test/150, corr1_true_test/150)
					print (corr2_test/150, corr2_true_test/150)
					print (corr3_test/150, corr3_true_test/150)
						
					

# execute for the three different types					
for type in range(NType):
	train_model(NArray[type],NameArray[type])
		

