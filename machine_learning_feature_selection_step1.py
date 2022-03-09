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
mode_array = ['basic']
#dyn_array = ['MD', 'LOG(FRES)']
#mode_array = ['basic', 'basic_phys_inh', 'basic_phys']

# read simulation data
KA2D_labels_full = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst/struct_filion_labels_type0.csv",nrows=NSTrain*NLOC)
KA2D_features_phys = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst/struct_filion_phys_type0.csv",nrows=NSTrain*NLOC)
KA2D_features_phys = KA2D_features_phys.iloc[: , :-1]
KA2D_features_phys=KA2D_features_phys.loc[:,KA2D_features_phys.columns.str.contains('INH')]
KA2D_features_phys=KA2D_features_phys.loc[:,~KA2D_features_phys.columns.str.contains('STD')]
KA2D_features_phys=KA2D_features_phys.loc[:,~KA2D_features_phys.columns.str.contains('PSI')]
print(KA2D_features_phys)
KA2D_features_thermal = pd.read_csv("../../../evaluatation/KA_model/eval_ml_T0.44_bapst/struct_filion_thermal_type0.csv",nrows=NSTrain*NLOC)

for m in mode_array:
	for d in dyn_array:
		for t in range(NT-3, NT):
			
			print ("Analyze mode ", m, " for observable ", d, " at time ", t)
		
			# choose labels
			KA2D_labels = KA2D_labels_full["{}{}".format(d, t)]

			if m=='basic' :
				KA2D_features = KA2D_features_thermal
				
			if m=='basic_phys_inh' or m=='basic_phys' :	
				KA2D_features = pd.concat([KA2D_features_phys,KA2D_features_thermal], axis=1)

			# split into learning and test set (for feature selection)
			sep=int(NSTrain/2)
			KA2D_features, KA2D_features_test = np.split(KA2D_features, [sep*NLOC])
			KA2D_labels, KA2D_labels_test = np.split(KA2D_labels, [sep*NLOC])

			step1=int((NSTrain-sep)/4)
			step2=int((NSTrain-sep)/2)
			step3=int((NSTrain-sep)/4*3)
			KA2D_features0, KA2D_features1, KA2D_features2, KA2D_features3 = np.split(KA2D_features, [step1*NLOC,step2*NLOC,step3*NLOC])
			KA2D_labels0, KA2D_labels1, KA2D_labels2, KA2D_labels3 = np.split(KA2D_labels, [step1*NLOC,step2*NLOC,step3*NLOC])
			
			# optimization loop
			corr0_max = -1.0
			names_list_max = []
			nstart=110
			nstep=10

			
			for x in range(100):
				rr0 = Ridge(alpha=lambda_ridge) 
				rr0.fit(KA2D_features0, KA2D_labels0)
				params0=rr0.coef_
				rr1 = Ridge(alpha=lambda_ridge) 
				rr1.fit(KA2D_features1, KA2D_labels1)
				params1=rr1.coef_
				rr2 = Ridge(alpha=lambda_ridge) 
				rr2.fit(KA2D_features2, KA2D_labels2)
				params2=rr2.coef_
				rr3 = Ridge(alpha=lambda_ridge) 
				rr3.fit(KA2D_features3, KA2D_labels3)
				params3=rr3.coef_

				params = np.array([params0, params1, params2, params3])
				
				mean=np.mean(params,axis=0)
				std=np.std(params,axis=0)
				score=std/abs(mean)
				score = (-score).argsort()
				indices1 = score[:nstart]
				remove_list = []
				if x == 0:
					indices= (abs(mean)).argsort()[:50]
					indices_total = np.concatenate((indices1, indices), axis=0)
					
					# remove features with smallest mean and large std to mean ratio
					remove_list = KA2D_features0.columns[indices_total]
					#print (remove_list)
					remove_list=[item for item in remove_list if 'DEN|EPOT' not in item]
					remove_list=[item for item in remove_list if 'DEN' not in item]	
					remove_list=[item for item in remove_list if 'EPOT' not in item]
				else :
					remove_list = KA2D_features0.columns[indices1]
					#print (remove_list)
					remove_list=[item for item in remove_list if 'DEN' not in item]	
					remove_list=[item for item in remove_list if 'EPOT' not in item]
					count = 1
					#print (len(remove_list))
					while len(remove_list) < nstart :
						remove_list.append(KA2D_features0.columns[score[nstart+count]])
						#print (remove_list)
						remove_list=[item for item in remove_list if 'DEN' not in item]	
						remove_list=[item for item in remove_list if 'EPOT' not in item]
						count +=1
						#print (len(remove_list))
					
	
				KA2D_features0.drop(remove_list, axis=1, inplace=True)
				KA2D_features1.drop(remove_list, axis=1, inplace=True)
				KA2D_features2.drop(remove_list, axis=1, inplace=True)
				KA2D_features3.drop(remove_list, axis=1, inplace=True)
				KA2D_features.drop(remove_list, axis=1, inplace=True)	
				KA2D_features_test.drop(remove_list, axis=1, inplace=True)	
				rr0 = Ridge(alpha=lambda_ridge) 
				rr0.fit(KA2D_features, KA2D_labels)
				
				#print(KA2D_features_test)
				y_out0=rr0.predict(KA2D_features_test)
				corr0, _ = pearsonr(y_out0, KA2D_labels_test)
				if corr0 > corr0_max:
					corr0_max = corr0
					# include names that were removed
					names_list_max = []
					for name in KA2D_features0.columns:
						names_list_max.append(name)
					
				length=KA2D_features0.shape[1]	
				print(x, length, corr0)
				if nstart > nstep :
					nstart -= nstep

				if (length < 100) :
					break;
					
			with open("features_{}{}{}.dat".format(m, d,t), 'w') as filehandle:
				for listitem in names_list_max:
					filehandle.write('%s\n' % listitem)


	

