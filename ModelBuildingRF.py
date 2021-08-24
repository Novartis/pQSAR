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


import os, sys
import math
import random
import scipy.stats
import pandas as pd
import numpy as np

from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

from sklearn.ensemble import RandomForestRegressor
import joblib
import CommonTools as Tools

CompID = 'compoundID'
smiles = 'smiles'
separator = '\t'
molecule = 'molecule'
MolCutoff = 50
fraction = 0.75
STDcutoff = 0.5
trees = 200
max_feature = 'sqrt'
seed = 42
#clusterCutoff = 0.43
assay_ver_id = 'AssayID'
PIC50 = 'PIC50'

#read optional arguments
if (len(sys.argv) > 5):
	fraction = float(sys.argv[5])
print("Training ratio: ", fraction)

if (len(sys.argv) > 6):
	MolCutoff = int(sys.argv[6])
print("Compound cutoff: ", MolCutoff)

if (len(sys.argv) > 7):
	CompID = sys.argv[7]
print("Compound ID: ", CompID)

if (len(sys.argv) > 8):
	smiles = sys.argv[8]
print("Smiles Field Name: ", smiles)

if (len(sys.argv) > 9):
	separator = sys.argv[9]
print("Field separator: Ignored")

if (len(sys.argv) > 10):
	trees = int(sys.argv[10])
print("Number of trees (estimator): ", trees)

if (len(sys.argv) > 11):
	njobs = int(sys.argv[11])



if os.path.isfile(sys.argv[3]):
	print('\tModel is already built')
	sys.exit()

data = pd.read_csv(sys.argv[1], header=0, index_col=0, sep=Tools.separator4fileName(sys.argv[1]))

col_names = list(data.columns)
col_names[0 : 3] = [assay_ver_id, PIC50, smiles]
data.columns = col_names
data = data[~data.index.duplicated()]
# For duplicates identified by smiles string, keep the most active one
data = data.sort_values(by=['PIC50']).drop_duplicates(subset=['smiles'], keep='last')

if (len(data) < MolCutoff):
	sys.exit()

try:
    PandasTools.AddMoleculeColumnToFrame(data,smiles,molecule)
except:
    print("Erroneous exception was raised and captured...")

#remove records with empty pIC50s
data = data.loc[data['PIC50'].notnull()]
#print(len(data))

#remove records with empty molecules
data = data.loc[data[molecule].notnull()]

data['FP'] = [Tools.computeFP(m) for m in data[molecule]]

#filter potentially failed fingerprint computations
data = data.loc[data['FP'].notnull()]
std = np.std(data[PIC50])

if (len(data) < MolCutoff) or (std < STDcutoff):
	sys.exit()

data.sort_index(axis=0, ascending=True, inplace=True)

if len(data) < 10_000:
	model = RandomForestRegressor(n_estimators=trees, max_features=max_feature, n_jobs=njobs, random_state=seed)
elif 10_000 <= len(data) < 50_000:
	model = RandomForestRegressor(n_estimators=trees, max_features=max_feature, min_samples_split=7, min_samples_leaf=4, n_jobs=njobs, random_state=seed)
else:
	model = RandomForestRegressor(n_estimators=trees, max_features=max_feature, min_samples_split=15, min_samples_leaf=8, n_jobs=njobs, random_state=seed)

#resolve wrapped fingerprints
x = [x.fp for x in data['FP']]
y = data[PIC50]
model.fit(x, y)

pIC50_min = np.min(data[PIC50])
pIC50_max = np.max(data[PIC50])

ActiveCount = len(data.loc[data[PIC50] >= 6])
ASSAY_VER_ID = data[assay_ver_id][0]

#save the model
joblib.dump(model, sys.argv[2])
print("Model File written")

#Getting model/dataset stats for random
random.seed(seed)

trainData = data.sample(frac=fraction, random_state=seed)
testData = data.drop(list(trainData.index), axis=0)

x = [x.fp for x in trainData['FP']]
y = trainData[PIC50]
model.fit(x, y)
testData['prediction'] = model.predict([x.fp for x in testData['FP']])

TestFraction_RF = -1

#split data into training and test set using clustered split

	
#Splited = Tools.ChallengedSplit(data, fraction, "Auto")
Splited = Tools.ChallengedSplit(data, smiles, fraction, "Auto", False, False)

print("Splits generated")
histogram = Tools.SimHistogram([AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in Splited.testSet[molecule]],[AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in Splited.trainSet[molecule]])
print("Histogram generated")
TestFraction_RF = round(float(len(Splited.trainSet))/float(len(data)),3)
print("Challenge set: training size: %s; test size: %s; fraction: %s" % (len(Splited.trainSet), len(Splited.testSet), TestFraction_RF))

data['Set'] = 'TRN'

data.loc[Splited.testSet.index, 'Set'] = 'TST'

x = [x.fp for x in Splited.trainSet['FP']]
y = Splited.trainSet[PIC50]
model.fit(x, y)
data['prediction'] = model.predict([x.fp for x in data['FP']])
if len(Splited.testSet) > 0:
	Splited.testSet['prediction'] = model.predict([x.fp for x in Splited.testSet['FP']])
	R2ext = round(math.pow(scipy.stats.pearsonr(Splited.testSet[PIC50], Splited.testSet['prediction'])[0], 2), 2)
else:
	R2ext = 'nan'


#Collecting the values
R2rand = round(math.pow(scipy.stats.pearsonr(testData[PIC50], testData['prediction'])[0], 2), 2)
R2fit = round(math.pow(scipy.stats.pearsonr(data[PIC50], data['prediction'])[0], 2), 2)
std = round(std, 2)

print("Assay: %s, R^2_RF(ext): %s, Molecules#: %s, R^2_RF(random): %s, STD (pIC50): %s, R^2_RF(fit): %s\n" %(ASSAY_VER_ID, R2ext, len(data), R2rand, std, R2fit))

#write stats and data to files
text_file = open(sys.argv[3], "w")
text_file.write("AssayID\tR^2_RF(ext)\tCount(molecules)\tR^2_RF(random)\tActives(>=6)\tpIC50_min\tpIC50_max\tSTD(pIC50)\tR^2_RF(fit)\tTrainFraction_RF\t")

#header
try:
	for i in range(1, len(histogram[1])):
		text_file.write("%s<->%s\t"%(histogram[1][i - 1], histogram[1][i]))
except NameError:
  print ("Histogram is not defined. The dataset is too big")
 
text_file.write("\n")

#values
text_file.write("%s\t%0.2f\t%i\t%0.2f\t%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t" %(ASSAY_VER_ID, R2ext, len(data), R2rand, ActiveCount, pIC50_min, pIC50_max, std, R2fit, TestFraction_RF))

try:
	for i in range(0, len(histogram[0])):
		text_file.write("%s\t"%(histogram[0][i]))
except NameError:
  print ("Omitting histogram")
text_file.write("\n")

text_file.close()

if (len(data) >= MolCutoff) and (std >= STDcutoff):
	#data['ASSAY_VER_ID'] = ASSAY_VER_ID
	data.drop([molecule, 'FP'], axis=1, inplace=True)
	try:
		data.drop(["Units", "ResultType"], axis=1, inplace=True)
	except:
		print("Fields Units and/or ResultType are not in the dataset!")

	data.to_csv(sys.argv[4], sep=Tools.separator4fileName(sys.argv[4]))
