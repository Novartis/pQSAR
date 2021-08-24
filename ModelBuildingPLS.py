"""
Copyright 2021 Novartis Institutes for BioMedical Research Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""


import os
import sys
import numpy as np
import pandas as pd
import random
from scipy.stats import pearsonr
from collections import namedtuple
from sklearn import cross_decomposition
import joblib
import CommonTools as Tools


#default values for some parameters
#os.environ['OPENBLAS_NUM_THREADS'] = '1'


def split_data(data, trn_idx):
	dat = data.set_index('compoundID')
	train = dat.loc[trn_idx]
	test = dat.drop(trn_idx, axis=0)
	train.reset_index(inplace=True)
	test.reset_index(inplace=True)
	Split = namedtuple('Split', 'testSet trainSet')
	return Split(testSet=test, trainSet=train)
	
	
# adjustable based on your dataset
MIN_NUMBER_COLUMNS = 10
##
##MIN_NUMBER_COLUMNS = 2
smiles = "smiles"
fraction = 0.75
assayName = sys.argv[4]


MaxCorrelation = 0.98
pIC50 = sys.argv[3]
fraction = float(sys.argv[8])
seed = 42

try:
	save_random = bool(sys.argv[9])
except:
	save_random = False

print("Training fraction: ", fraction)

if os.path.isfile(sys.argv[5]):
	print('\tModel is already there ...')
	sys.exit()


if os.path.isdir(sys.argv[1]):
	#Read file, but sometimes a header is missing
	data = pd.read_csv(sys.argv[2], sep=Tools.separator4fileName(sys.argv[2]), index_col=0, header=0, error_bad_lines=False)

	#name of the assay column
	col_names = list(data.columns)
	col_names[0 : 4] = ['AssayID', pIC50, smiles, 'Set']
	data = data[col_names[0 : 4]]
	data.drop('AssayID', axis=1, inplace=True)
	NVS_IDs = list(data.index)

	p = Tools.QueryCmpd(sys.argv[1])
	p_col = p.columns()

	df_cmpd = pd.DataFrame(0.0, index=NVS_IDs, columns=p_col, dtype=np.float32)

	for cmpd in NVS_IDs:
		#print(cmpd, thr)
		try:
			name, pqsar_vec = p.get(cmpd)
			df_cmpd.loc[cmpd] = [np.float32(s) for s in pqsar_vec]
		except:
			df_cmpd.drop(cmpd, axis=0, inplace=True)
			data.drop(cmpd, axis=0, inplace=True)
			print('Warning! {} not found'.format(cmpd))

	data = df_cmpd.merge(data, how='left', right_index=True, left_index=True)
	data.index.name = 'compoundID'
	trn_idx = list(data.loc[data['Set'] == 'TRN'].index)
	trn_tst_list = list(data['Set'])
	data.drop('Set', axis=1, inplace=True)
	data.reset_index(inplace=True)
	
	del df_cmpd
	del p
else:
	data = pd.read_csv(sys.argv[2], sep=Tools.separator4fileName(sys.argv[2]), index_col=0, header=0, error_bad_lines=False)
	data.index.name = 'compoundID'
	data.sort_index(axis=0, ascending=True, inplace=True)
	
#should always be True!
is2removeColumn = True

#clean up in case of irregular rows
print ("Before dataset clean-up: ", len(data))

float_cols = [c for c in data if data[c].dtype == "float64"]
#print(num_columns)
print("Total assay columns: ", len(float_cols))

#drop columns with missing numerical values
### using downcast would convert datatype from float64 to float32. save half of the memory usage.
#data[num_columns] = data[num_columns].apply(pd.to_numeric, downcast='float', errors='coerce')

#print("Here")

data.dropna(axis=0, subset=float_cols, inplace=True)
data[float_cols] = data[float_cols].apply(pd.to_numeric, downcast='float', errors='coerce')
#remove rows with IDs that are only one character long
#data = data[data.iloc[:,0].str.len() > 1]

if assayName in data.columns:
	data.drop(assayName, axis=1, inplace=True)
	
print("After dataset clean-up: ", len(data))

if os.path.isdir(sys.argv[1]):
	Splited = split_data(data, trn_idx)

elif sys.argv[1] == 'custom':
	trn_idx = list(data.loc[data['Set'] == 'TRN'].index)
	trn_tst_list = list(data['Set'])
	data.index.name = 'compoundID'
	data.reset_index(inplace=True)
	data['smiles'] = 'NA'
	data.drop('Set', axis=1, inplace=True)
	Splited = split_data(data, trn_idx)
	data.reset_index(inplace=True)
	
else:
	data.reset_index(inplace=True)
	Splited = Tools.ChallengedSplit(data, smiles, fraction)
	trn_idx = list(Splited.trainSet.iloc[:, 0])
	cmpd_list = data[['compoundID']].copy()
	cmpd_list.set_index('compoundID', inplace=True)
	cmpd_list['Set'] = 'TST'
	cmpd_list.loc[trn_idx] = 'TRN'
	trn_tst_list = list(cmpd_list['Set'])
	


std = np.std(data[pIC50])
r2fit = None
q2orig = None
nv = None
pred_fit = None
r2ex_bar = -1
dscr = []
THRESHOLD = 0
mae = np.inf
best_model = ''
###
#split data into training and test set using clustered split
if len(data) < 30_000:
	MinCorrelation = [0.001, 0.05, 0.20]
	### cutoff 0.001 to remove nearly half variables for extremely large datasets without harming model performance
else:
	MinCorrelation = [0.20]

print("{}: Challenged split created. {} in the training set and {} molecules in the test set".format(assayName, len(Splited.trainSet), len(Splited.testSet)))

###calculate the correlation based on the training Y data.
Corr = Tools.data_corr(Splited.trainSet, assayName, pIC50)
Corr_set = list(set(Corr.isnull().all()))

#if not len(Corr_set) == 1 and Corr_set[0] == True:
if not (len(Corr_set) == 1 and Corr_set[0] == True):
	for i, thr in enumerate(MinCorrelation):
		profile = [v for v in Corr.columns if float(Corr[v]) >= thr and float(Corr[v]) <= MaxCorrelation]
		profile = list(set(profile).intersection(set(data.columns)))

		print ("%i numeric columns in profile"%(len(profile)))
		
		if i == 0:
			while len(profile) < MIN_NUMBER_COLUMNS:
				thr = 3 * thr / 4
				profile = [v for v in Corr.columns if float(Corr[v]) >= thr and float(Corr[v]) <= MaxCorrelation]
				profile = list(set(profile).intersection(set(data.columns)))

		if len(profile) < MIN_NUMBER_COLUMNS:
			break

		profile.extend([pIC50, smiles])
		profile.insert(0, 'compoundID')
		trainSet = Splited.trainSet[profile]
		testSet = Splited.testSet[profile]

		#build the model based on that training split
		try:
			PLSchallenged = Tools.buildPLS(trainSet, pIC50, smiles)
		except:
			# in case assays with only a few activies.
			PLSdataset = Tools.extractPLSparts(trainSet, pIC50, smiles)
			# randomly chosen
			nv = 8
			model = cross_decomposition.PLSRegression(n_components=nv, scale=False).fit(PLSdataset.independent, PLSdataset.dependent)
			q2orig = np.nan
			PLS = namedtuple('PLS', 'model q2orig nv')
			PLSchallenged = PLS(model=model, q2orig=q2orig, nv=nv)



		print("Challenged model built from %i molecules" % (len(trainSet)))

		#extract independent variables from the test part of the split for doing predictions on them
		PLSdataset = Tools.extractPLSparts(testSet, pIC50, smiles)
		print("Independent variables extracted from challenged model for predictions")

		#calculate predictions on the test set using the model built on the test set
		pred_challenged = PLSchallenged.model.predict(PLSdataset.independent)
		print("Test set predictions for challenged model completed")

		#Calculate r2ext - predicted on the partial model using the test set
		r2ex = pearsonr(PLSdataset.dependent, pred_challenged.transpose()[0])
		r2ex = round(r2ex[0]**2,3)
	
		if np.isnan(r2ex):
			# mean absolute error (mae) as the performance indicator if R2ext is not available
			pred_mae = sum(abs(PLSdataset.dependent - pred_challenged.transpose()[0])) / len(PLSdataset.dependent)
			if mae >= pred_mae:
				mae = pred_mae
				dscr = profile
				r2ex_bar = r2ex
				THRESHOLD = thr
				best_model = PLSchallenged
		else:	
			if r2ex_bar <= r2ex:
				r2ex_bar = r2ex
				dscr = profile
				THRESHOLD = thr
				best_model = PLSchallenged

	###
	# to save memory. please remove these lines if not elegant
	del Splited
	try:
		del trainSet
	except:
		print("trainSet not defined!") 
	try:
		del testSet
	except:
		print("testSet not defined!") 
	try:
		del PLSchallenged
	except:
		print("PLSchallenged not defined!") 
	try:
		del PLSdataset
	except:
		print("PLSdataset not defined!") 
	
	###
	r2ex = r2ex_bar

else:
	dscr = list(Corr.columns)
	dscr.extend([pIC50, smiles])
	dscr.insert(0, 'compoundID')
	r2ex = np.nan


#if model file exists, do not redo the model
if not os.path.isfile(sys.argv[5]):
	
	### reduce profile to the best
	data = data[dscr]

	#build PLS model on all data
	try:
		PLSall = Tools.buildPLS(data, pIC50, smiles)
	except:
		PLSdataset = Tools.extractPLSparts(data, pIC50, smiles)
		# randomly chosen
		nv = 8
		model = cross_decomposition.PLSRegression(n_components=nv, scale=False).fit(PLSdataset.independent, PLSdataset.dependent)
		q2orig = np.nan
		PLS = namedtuple('PLS', 'model q2orig nv')
		PLSall = PLS(model=model, q2orig=q2orig, nv=nv)

	#extract independent variables from the large set to calculate predictions and compare to experimental data
	PLSall_split = Tools.extractPLSparts(data, pIC50, smiles)

	#calculate predictions on all molecules using the model built from all molecules
	if best_model != '':
		pred_fit = best_model.model.predict(PLSall_split.independent)
	else:
		pred_fit = PLSall.model.predict(PLSall_split.independent)
	
	if save_random == True:
		random.seed(seed)
		train_random = data.sample(frac=fraction, random_state=seed).copy()
		test_random = data.drop(train_random.index, axis=0)
		PLS_random = Tools.buildPLS(train_random, pIC50, smiles)
		test_random = Tools.extractPLSparts(test_random, pIC50, smiles)
		random_test = PLS_random.model.predict(test_random.independent).transpose()[0]
		r2_random = round(pearsonr(test_random.dependent, random_test)[0] ** 2, 3)
		random_fit = PLS_random.model.predict(PLSall_split.independent).transpose()[0]
		
		trn_random_list = list(train_random.iloc[:, 0])

	#Calculate r2ext - predicted on the partial model using the test set
	r2fit = pearsonr(PLSall_split.dependent,pred_fit.transpose()[0])
	r2fit = round(r2fit[0]**2,3)
	q2orig = PLSall.q2orig
	nv = PLSall.nv

	#save it
	
	joblib.dump((PLSall.model, r2fit, PLSall.q2orig, PLSall.nv, dscr, THRESHOLD), sys.argv[5])
	print("Wrote model file and model parameters to the same file: %s" % sys.argv[5])
	
else:
	saved = joblib.load(open(sys.argv[5], 'rb'))
	r2fit = saved[1]
	q2orig = saved[2]
	nv = saved[3]
	print("Model file %s exists and will not be overridden" % sys.argv[5])

if q2orig is not None:
	#print ("Assay: %s; r2ex = %0.2f; q2orig = %0.2f; nv = %i; r2fit = %0.2f; std(pIC50) =  %0.2f" % (assayName, r2ex, PLSall.q2orig, PLSall.nv, r2fit, std))
	print ("Assay: %s; r2ex=%0.2f; q2orig=%0.2f; nv=%i; r2fit=%0.2f; std(pIC50)= %0.2f" % (assayName, r2ex, q2orig, nv, r2fit, std))
	
	#Save summary
	#summary = pd.DataFrame({"AssayID":assayName, "R2ext":r2ex, "NV":PLSall.nv, "Q2orig":PLSall.q2orig, "R^2fit":r2fit, "std(pIC50)":std, "#compounds":len(data)},index = [0])
	if save_random == True:
		summary = pd.DataFrame({"AssayID":assayName, "R2ext":r2ex, "NV":nv, "Q2orig":q2orig, "R^2fit":r2fit, "R2random": r2_random,"std(pIC50)":std, "Compounds#":len(data), "Threshold":THRESHOLD, "Columns#":(len(dscr) - 3)},index = [0])
		
	else:
		summary = pd.DataFrame({"AssayID":assayName, "R2ext":r2ex, "NV":nv, "Q2orig":q2orig, "R^2fit":r2fit, "std(pIC50)":std, "Compounds#":len(data), "Threshold":THRESHOLD, "Columns#":(len(dscr) - 3)},index = [0])
		
		#due to a bug, it sorts columns in an alphabetical order
	summary.to_csv(sys.argv[6], index=False, sep=Tools.separator4fileName(sys.argv[6]), float_format="%.2f")
	
else: 
	print ("Assay: %s; r2ex=%0.2f; q2orig=%0.2f; nv=%i; r2fit=%0.2f; std(pIC50)=%0.2f" % (assayName, r2ex, -1, -1, r2fit, std))

if pred_fit is not None:
	#Save predicted vs experimental results for experimental compounds
	predVexp = pd.concat([data.iloc[:,0], data.loc[:, pIC50], pd.DataFrame(pred_fit,columns=["pIC50_PLS"], index=data.index)], axis=1 )
	predVexp.columns = ['compoundID','pIC50', 'pIC50_PLS']
	predVexp.set_index('compoundID', inplace=True)
	predVexp['AssayID'] = assayName
	predVexp['Set'] = trn_tst_list
	if save_random == True:
		predVexp['pIC50_random_pred'] = random_fit
		predVexp['randomSet'] = 'TST'
		predVexp.loc[trn_random_list, 'randomSet'] = 'TRN'

	predVexp.to_csv(sys.argv[7], index=True, header=True, sep=Tools.separator4fileName(sys.argv[7]), float_format="%.3f")
	
