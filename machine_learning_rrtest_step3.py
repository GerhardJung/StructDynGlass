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
NSTtest=70

lambda_ridge=0.2

dyn_array = ['MD']
mode_array = ['basic']
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
			KA2D_labels_full = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst/struct_filion_labels_type0.csv",nrows=(NSTrain+NSTtest)*NLOC, skiprows=[i for i in range(1,NSTrain*NLOC+1)])
			KA2D_features_phys = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst/struct_filion_phys_type0.csv",nrows=(NSTrain+NSTtest)*NLOC, usecols=names_list_phys, skiprows=[i for i in range(1,NSTrain*NLOC+1)])
			#print(KA2D_features_phys)
			KA2D_features_thermal = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst/struct_filion_thermal_type0.csv",nrows=(NSTrain+NSTtest)*NLOC, usecols=names_list_cg, skiprows=[i for i in range(1,NSTrain*NLOC+1)])
			#print(KA2D_features_thermal)

		
			# choose labels
			KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]

			if m=='basic' :
				KA2D_features = KA2D_features_thermal
				
			if m=='basic_phys_inh' or m=='basic_phys' :	
				KA2D_features = pd.concat([KA2D_features_phys,KA2D_features_thermal], axis=1)
					
			# read learned weights
			filehandle = open("rrparams_{}{}{}.dat".format(m, d, t), "r")
			params = filehandle.read()
			params_list = params.split("\n")
			del params_list[-1]
			inter=float(params_list[-1])
			del params_list[-1]
			params_list_np = np.zeros(len(params_list), dtype=np.float32)
			for x in range(len(params_list)):
				params_list_np[x] = float(params_list[x])
			#print (params_list_np)
			
			#y_pred = np.random.uniform(low=-0.01, high=0.01, size=NSTtest*NLOC)
			y_pred = np.zeros(NSTtest*NLOC)
			#print(y_pred)
			#print(KA2D_features.values)
			#print(y_pred.shape,params_list_np.shape, KA2D_features.values.shape)
			
			y_pred = np.matmul(KA2D_features.values, params_list_np)
			y_pred += float(inter)
			
			corr0, _ = pearsonr(y_pred, KA2D_labels)

			print(corr0)


	

