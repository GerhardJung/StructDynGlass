from signal import SIGUSR1
import pandas as pd
import numpy as np

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
import statistics

import tensorflow_probability as tfp
import tensorflow as tf
from tensorflow.python.ops import math_ops
from tensorflow.keras import layers
from keras.layers import Dropout
from tensorflow.keras import optimizers
from tensorflow.keras import initializers
from tensorflow.keras import regularizers
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

N=4096
NType=2
NLOC=3277
#N1=819
NT=10
NSTrain=70
NSTest=70

costfct_lambda=0.0

dyn_array = ['MD']
mode_array = ['basic_phys_inh']
#dyn_array = ['MD', 'LOG(FRES)']
#mode_array = ['basic', 'basic_phys_inh', 'basic_phys']

for m in mode_array:
	for d in dyn_array:
		for t in range(NT-1, NT):
			
			print ("Analyze mode ", m, " for observable ", d, " at time ", t)

			# first read in labels
			filehandle = open("features_{}{}{}.dat".format(m, d,t), "r")
			names = filehandle.read()
			names_list = names.split("\n")
			del names_list[-1]
			#print(names_list)
			names_list_phys= []
			names_list_cg= []
			for x in names_list:
				if x[:2] == "CG":
					names_list_cg.append(x)
				else :
					names_list_phys.append(x)

			# read simulation data
			KA2D_labels_full = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst2/struct_filion_labels_type0.csv",nrows=(NSTrain+NSTest)*NLOC)
			KA2D_features_phys = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst2/struct_filion_phys_type0.csv",nrows=(NSTrain+NSTest)*NLOC, usecols=names_list_phys)
			#print(KA2D_features_phys)
			KA2D_features_thermal = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst2/struct_filion_thermal_type0.csv",nrows=(NSTrain+NSTest)*NLOC, usecols=names_list_cg)
			#print(KA2D_features_thermal)

		
			# choose labels
			KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]

			if m=='basic' :
				KA2D_features = KA2D_features_thermal
				
			if m=='basic_phys_inh' or m=='basic_phys' :	
				KA2D_features = pd.concat([KA2D_features_phys,KA2D_features_thermal], axis=1)
				
			length=KA2D_features.shape[1]
			print(length)
			
			KA2D_features, KA2D_features_test = np.split(KA2D_features, [NSTrain*NLOC])
			KA2D_labels, KA2D_labels_test = np.split(KA2D_labels, [NSTrain*NLOC])
					
			# create ML model
			Nbottle = 5
			Ninner = 20
			KA2D_model = tf.keras.Sequential()
			ReLUlayer = tf.keras.layers.ELU(alpha=5.0)
			KA2D_model.add(
			  Dense(Nbottle, 
			  activation=ReLUlayer, 
			  input_shape=(length,),
			  #kernel_initializer=initializers.RandomNormal(stddev=0.2),
			  bias_initializer=initializers.RandomNormal(stddev=0.0),		  
        		  kernel_regularizer=regularizers.l2(costfct_lambda)
        		  )
			) 
			KA2D_model.add(
			  Dense(Ninner,activation=ReLUlayer,
			  #kernel_initializer=initializers.RandomNormal(stddev=0.3),
			  bias_initializer=initializers.RandomNormal(stddev=0.0),		  
        		  kernel_regularizer=regularizers.l2(costfct_lambda)
			  )
			) 
			KA2D_model.add(
			  Dense(Ninner,activation=ReLUlayer,
			  #kernel_initializer=initializers.RandomNormal(stddev=0.3),
			  bias_initializer=initializers.RandomNormal(stddev=0.0),		  
        		  kernel_regularizer=regularizers.l2(costfct_lambda)
			  )
			) 
			KA2D_model.add(
			  Dense(1,activation='linear',		  
        		  kernel_regularizer=regularizers.l2(costfct_lambda)
			  )
			) 
			print(KA2D_model.summary())

			KA2D_model.compile(loss='mse',optimizer = optimizers.Adam(0.002))
			#print(KA2D_model.layers[0].get_weights()    )
			for x in range(25):
				history = KA2D_model.fit(
				    KA2D_features,
				    KA2D_labels,
				    batch_size=500,
				    epochs=10,
				   )

				
				y_out=KA2D_model.predict_on_batch(KA2D_features_test)[:,0]
				corr, _ = pearsonr(y_out, KA2D_labels_test)
				print (x, corr)
	

