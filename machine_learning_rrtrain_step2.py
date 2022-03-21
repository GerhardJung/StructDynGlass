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

N=4096
NType=2
NLOC=3277
#N1=819
NT=10
NSTrain=70

lambda_ridge=0.2

dyn_array = ['MD']
mode_array = ['basic_phys_inh']
#dyn_array = ['MD', 'LOG(FRES)']
#mode_array = ['basic', 'basic_phys_inh', 'basic_phys']

for m in mode_array:
	for d in dyn_array:
		for t in range(NT-3, NT):
			
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
			KA2D_labels_full = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst2/struct_filion_labels_type0.csv",nrows=NSTrain*NLOC)
			KA2D_features_phys = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst2/struct_filion_phys_type0.csv",nrows=NSTrain*NLOC, usecols=names_list_phys)
			#print(KA2D_features_phys)
			KA2D_features_thermal = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst2/struct_filion_thermal_type0.csv",nrows=NSTrain*NLOC, usecols=names_list_cg)
			#print(KA2D_features_thermal)

		
			# choose labels
			KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]

			if m=='basic' :
				KA2D_features = KA2D_features_thermal
				
			if m=='basic_phys_inh' or m=='basic_phys' :	
				KA2D_features = pd.concat([KA2D_features_phys,KA2D_features_thermal], axis=1)
				
			length=KA2D_features.shape[1]
			print(length)
					
			rr0 = Ridge(alpha=lambda_ridge) 
			rr0.fit(KA2D_features, KA2D_labels)
			
			params=rr0.coef_	
			inter=rr0.intercept_

			with open("rrparams_{}{}{}.dat".format(m, d,t), 'w') as filehandle:
				for listitem in params:
					filehandle.write('%f\n' % listitem)
				filehandle.write('%f\n' % inter)


	

