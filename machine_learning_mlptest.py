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
import matplotlib
matplotlib.rcParams['figure.dpi']=300 # highres display

from scipy.stats import pearsonr
from scipy.stats import spearmanr
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

Temp=0.251
NT=42
NSTrain=150
NSTest=150
Tstart=25
Tstop=NT

dyn_array = ['BB']
mode_array = ['basic_phys_inh']
#dyn_array = ['MD', 'LOG(FRES)']
#mode_array = ['basic', 'basic_phys_inh', 'basic_phys']

# data to save particle predictions
#NTsave=int(Tstop-Tstart)
#Nmodes=len(dyn_array)*len(mode_array)
#NCG=3
#data_pred_type0 = np.empty((NTsave*Nmodes*NCG,NSTest*NLOC), dtype=float)

# read time data
time_data = np.empty((NT+30), dtype=float)
iso_data = np.empty((NT+30), dtype=float)
f = open("../../../evaluatation/KA2D_model/eval_ml_T" + str(Temp) +"/isoconf_predictability_BB.dat", "r")
count = 1
for x in f:
	time_data[count]=x.split()[0]
	iso_data[count]=x.split()[1]
	count+=1
	
def test_model(NLOC,Type):

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
	PATH='./results_mlp'
	isExist = os.path.exists(PATH)
	if not isExist:
	  os.makedirs(PATH)

	#count=0
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
				sourceFile = open(PATH + "/pearson_{}{}_{}_{}.dat".format(m, d, Type,cg), 'w')
				tcount=0
				for t in range(Tstart, Tstop):
					
					print ("Analyze ", Type, " mode ", m, " for observable ", d, " cg length ", cg," at time ", t)

										# choose labels
					KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]
					print(KA2D_labels)

					KA2D_features = KA2D_features_phys
					KA2D_features, KA2D_features_test = np.split(KA2D_features, [NSTrain*NLOC])
					KA2D_labels, KA2D_labels_test = np.split(KA2D_labels, [NSTrain*NLOC])
						
					length=KA2D_features.shape[1]
					print(length)
					
					# shuffle arrays
					#KA2D_features = KA2D_features.sample(frac=1,random_state=12345)
					#KA2D_labels = KA2D_labels.sample(frac=1,random_state=12345)
					#select labels
					NumSelect = int((NLOC*NSTrain)/5)
					smallest_index = KA2D_labels.nsmallest(NumSelect).index
					largest_index = KA2D_labels.nlargest(NumSelect).index
					select_index=smallest_index.union(largest_index)
					#print (largest_index)
					#KA2D_labels = KA2D_labels.iloc[select_index-2]
					print(KA2D_labels)
					#KA2D_features = KA2D_features.iloc[select_index-2]
					
					# load ML model
					def loss_carrier(factor1, factor2, factor3,NLOC_loc):
						def loss(y_true, y_pred):
						    	# loss1, squared distances between labels
							squared_difference = tf.abs(y_true[:,0] - y_pred[:,0])
							
							# calculate cg labels
							Mean = tf.reduce_mean(y_pred[:,0], axis=-1)
							Fluct = tf.math.multiply(tf.cast(y_pred[:,0]-Mean, tf.float32),tf.cast(y_pred[:,0]-Mean, tf.float32))
							Fluct = tf.reduce_sum(Fluct)
							res1=tf.linalg.matvec(y_true[:,7:NLOC+7],y_pred[:,0]-Mean)/y_true[:,4]
							res2=tf.linalg.matvec(y_true[:,7+NLOC:2*NLOC+7],y_pred[:,0]-Mean)/y_true[:,5]
							res3=tf.linalg.matvec(y_true[:,7+2*NLOC:3*NLOC+7],y_pred[:,0]-Mean)/y_true[:,6]
							
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
							
							#tf.print("\n res:", res1, res2, res3, output_stream=sys.stdout)
							#tf.print("\n true:", y_true[0,1], y_true[0,2], y_true[0,3], output_stream=sys.stdout)
							
							# loss2, squared distance between cg labels
							squared_difference1 = tf.abs(y_true[0,1] - res1)
							squared_difference2 = tf.abs(y_true[0,2] - res2)
							squared_difference3 = tf.abs(y_true[0,3] - res3)
							
							#tf.print("\n diff:", squared_difference1, squared_difference2, squared_difference3, output_stream=sys.stdout)

							#res1_loc=tf.linalg.matvec(wlist1[0],y_true[:,0])
							#tf.print("\n y_true_pre:", y_true[:,3], output_stream=sys.stdout)
							#tf.print("\n y_true_calc:", res3, output_stream=sys.stdout)
							#tf.print("\n y_true:", res1_loc, output_stream=sys.stdout)
					#		res1 = tf.math.multiply(res1,y_true[:,0]-mean)
					#		res1 = tf.reduce_sum(res1)
					#		res1 /= norm1[0]
							#tf.print("\n y_true:", tf.reduce_mean(y_true), output_stream=sys.stdout)
							#tf.print("\n val:", res1/norm*NLOC, output_stream=sys.stdout)
							
							loss = tf.reduce_mean(squared_difference, axis=-1) + factor1*squared_difference1 + factor2*squared_difference2 + factor3*squared_difference3
							#loss = tf.reduce_mean(squared_difference, axis=-1) 
							
							return loss
						return loss
					def pearson_loss(y_true, y_pred):
						loss=tfp.stats.covariance((y_true[:,0])[:,None],(y_pred[:,0])[:,None])
						norm1= tf.math.sqrt(tfp.stats.covariance((y_true[:,0])[:,None]))
						norm2= tf.math.sqrt(tfp.stats.covariance((y_pred[:,0])[:,None]))
						return loss/(norm1*norm2)
					KA2D_model = load_model("./models_mlp/mlpparams_{}{}{}_{}_{}.dat".format(m, d,t, Type,cg), custom_objects={'loss': loss_carrier(0.0,0.0,0.1,NLOC), 'pearson_loss': pearson_loss})
					   
					y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
					corr, _ = pearsonr(y_out, KA2D_labels_test)
					coef, _ = spearmanr(y_out, KA2D_labels_test)
					
					# save data
					#if t%4 == 0:
					#	data[count*NTsave*Nmodes+tcount,:] = y_out
					#	tcount+=1
					
					# calc chi_4
					var = mean = 0
					var_single = mean_single = 0
					for s in range(NSTest):
						qtot_true = 0
						labels = KA2D_labels[s*NLOC:(s+1)*NLOC]
						for x in labels :
							if d=='BB' :
								qtot_true+=x
								mean_single += x
								var_single += x*x
							else :
								if x > 1.0 :
									qtot_true+=1
									mean_single += 1
									var_single += 1
						#print(s,qtot_pred)
						mean += qtot_true
						var += qtot_true*qtot_true
					mean/=NSTest
					var/=NSTest
					mean_single/=NSTest*NLOC
					var_single/=NSTest*NLOC
					#print((var-mean*mean)/NLOC,var_single-mean_single*mean_single)
					chi4_true=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
					var = mean = 0
					var_single = mean_single = 0
					for s in range(NSTest):
						qtot_pred = 0
						labels = y_out[s*NLOC:(s+1)*NLOC]
						for x in labels :
							if d=='BB' :
								qtot_pred+=x
								mean_single += x
								var_single += x*x
							else :
								if x > 1.0 :
									qtot_pred+=1
									mean_single += 1
									var_single += 1
						#print(s,qtot_pred)
						mean += qtot_pred
						var += qtot_pred*qtot_pred
					mean/=NSTest
					var/=NSTest
					mean_single/=NSTest*NLOC
					var_single/=NSTest*NLOC
					chi4_pred=(var-mean*mean)/NLOC/(var_single-mean_single*mean_single)
					
					print (time_data[t], corr,coef, chi4_true, chi4_pred)
					print (time_data[t], corr,coef, chi4_true, chi4_pred, file = sourceFile)
					
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
						bins = np.linspace(0.75, 1, 100)
					
					#plt.hist([y_outrest,y_outlow,y_outhigh], bins,  alpha=0.5, label=['model rest','model low','model high'], stacked=True)
					plt.figure()
					plt.hist(KA2D_labels_test, bins,  alpha=0.5, label='true')
					plt.hist([y_outm2,y_outm1,y_out0,y_outp1,y_outp2], bins,  alpha=0.5, label=['model <0.2','model <0.4','model <0.6','model <0.8','model <1.0'], histtype='step', stacked=True, fill=False)
					plt.legend(loc='upper left')
					plt.title("{} {} t={}".format(m, d, time_data[t]))
					plt.xlabel("{}".format(d))
					plt.ylabel("PDF")
					plt.savefig(PATH + "/histogram_{}{}{}_{}_{}.png".format(m, d, t, Type,cg))
					
					plt.figure()
					plt.scatter(KA2D_labels_test, y_out,s=1,  c="g", alpha=0.5)
					plt.title("{} {} t={}".format(m, d, time_data[t]))
					plt.xlabel("true")
					plt.ylabel("pred")
					plt.savefig(PATH + "/scatter_{}{}{}_{}_{}.png".format(m, d, t, Type,cg))
					
					size=30.0
					iso_true=np.mean(KA2D_labels_test)
					iso_pred=np.mean(y_out)
					if Type == 'type0' :
						plt.figure()
						plt.title("{} t={}\n iso_true={:.2f} iso_pred={:.2f}".format(d, time_data[t], iso_true, iso_pred))
						plt.xlabel("X")
						plt.ylabel("Y")
						#print(plt.get_fignums())
						plt.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=size, cmap='coolwarm')
						
					if Type == 'type1' :
						ind = int(4.0 - cg)*4*(Tstop-Tstart) + 4*(t-Tstart) + 3
						#print(ind)
						plt.figure(ind)
						plt.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=0.94*size, cmap='coolwarm')

					if Type == 'type2' :
						ind = int(4.0 - cg)*4*(Tstop-Tstart) + 4*(t-Tstart) + 3
						#print(ind)
						plt.figure(ind)
						plt.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=y_out[-NLOC:],s=0.88*size, cmap='coolwarm')
						plt.colorbar()
						plt.savefig(PATH + "/conf_pred_{}{}{}_{}.png".format(m, d, t, cg))
						
					if Type == 'type0' :
						plt.figure()
						plt.title("{} t={}\n iso_true={:.2f}".format(d, time_data[t], iso_true))
						plt.xlabel("X")
						plt.ylabel("Y")
						#print(plt.get_fignums())
						plt.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=KA2D_labels_test[-NLOC:],s=size, cmap='coolwarm')
						
					if Type == 'type1' :
						ind = int(4.0 - cg)*4*(Tstop-Tstart) + 4*(t-Tstart) + 4
						#print(ind)
						plt.figure(ind)
						plt.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=KA2D_labels_test[-NLOC:],s=0.94*size, cmap='coolwarm')

					if Type == 'type2' :
						ind = int(4.0 - cg)*4*(Tstop-Tstart) + 4*(t-Tstart) + 4
						#print(ind)
						plt.figure(ind)
						plt.scatter(XPOS[-NLOC:],YPOS[-NLOC:], c=KA2D_labels_test[-NLOC:],s=0.88*size, cmap='coolwarm')
						plt.colorbar()
						plt.savefig(PATH + "/conf_true_{}{}{}.png".format(m, d, t))


				sourceFile.close()
	#			count += 1
			
	# print predicted configurations
	#config_file = open(PATH +"/pred_propensity_" + Type + ".dat", "w")
	#for m in mode_array:
	#	for d in dyn_array:
	#		for t in range(Tstart,Tstop):
	#			if t%4 == 0:
	#				config_file.write("{}{}{}_PRED ".format(m, d,t))
	#				config_file.write("{}{}{}_TRUE ".format(m, d,t))
	#config_file.write("\n") 
	#for p in range(NSTest*NLOC):
	#	count=0		
	#	for m in mode_array:
	#		for d in dyn_array:
	#			tcount=0
	#			for t in range(Tstart,Tstop):
	#				if t%4 == 0:	
	#					#print(t,count,count*NTsave*Nmodes+tcount)
	#					config_file.write("%f " % data[count*NTsave*Nmodes+tcount,p])
	#					config_file.write("%f " % KA2D_labels_full["{}{}".format(d, t)][p+2])
	#
	#					tcount += 1
	#			count += 1
	#	config_file.write("\n")
					
# execute for the three different types					
for type in range(NType):
	test_model(NArray[type],NameArray[type])
		

