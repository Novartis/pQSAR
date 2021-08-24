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
from scipy.spatial.distance import cdist, pdist
from scipy.cluster import hierarchy
import numpy as np
import joblib
from sklearn import cross_decomposition
from sklearn.model_selection import KFold
import pandas as pd

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

from collections import namedtuple


def data_corr(data, assayName, pIC50):
	mat = data.copy()
	expr = mat[[pIC50]]
	mat.drop([pIC50, mat.columns[0], 'smiles'], axis=1, inplace=True)
	stats = (1 - cdist(expr.T, mat.T, metric='correlation')) ** 2
	stats = pd.DataFrame(stats, index=[assayName], columns=mat.columns)
	return(stats)

def separator4fileName(File):
	return ',' if File.endswith('.csv') else '\t'

def extractPLSparts	(data, pIC50, smiles):
	#get everything from the column #1
	X = data.iloc[:, 1:]

	#what to do with the original column that is to be fitted? Normally it should be removed
	X.pop(pIC50)

	#remove smiles column
	if (smiles != None):
		X.pop(smiles)

	#format for PLS
	X = X.to_numpy()

	y = None
	try:
		#get column that is to be fitted (dependent variable)
		y = data.loc[:, pIC50]
		#also format it
		y = y.to_numpy()
	except:
		print("y will be empty, because there dependent column")

	PLSdataset = namedtuple('PLSdataset', 'dependent independent')
	return PLSdataset(dependent = y, independent = X)


def buildPLS(data, pIC50, smiles):

	PLSdataset = extractPLSparts(data, pIC50, smiles)
	y = PLSdataset.dependent
	X = PLSdataset.independent
	seed = 42
	q2 = -1
	nv = 1
	q2orig = -2
	cv = 5
	
	pc = X.shape[1] if X.shape[1] < 25 else 25 
	
	#Cross Validation
	kf = KFold(n_splits=cv, shuffle=True, random_state=seed)

	for n_comp in range(1, pc + 1):
		r2 = 0
		for train,test in kf.split(X):
			model = cross_decomposition.PLSRegression(n_components=n_comp, scale=False).fit(X[train], y[train])
			r2 += model.score(X[test], y[test])
		q2o = r2 / cv
		if q2o - (0.002*n_comp - 0.002) > q2:
			# penalty on the number of latent variables
			q2 = q2o - (0.002 * n_comp - 0.002)
			q2orig = q2o
			nv = n_comp

	#Model Building
	model = cross_decomposition.PLSRegression(n_components=nv,scale=False).fit(X,y)

	PLS = namedtuple('PLS', 'model q2orig nv')

	return PLS(model=model, q2orig=q2orig, nv=nv)
	#return {'model':model, 'q2':q2 ,'q2orig':q2orig, 'nv': nv }

class FP:
  def __init__(self, fp):
        self.fp = fp
  def __str__(self):
      return self.fp.__str__()

def computeFP(x):
	#compute depth-2 morgan fingerprint hashed to 1024 bits
	try:
		fp = AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024)
		res = np.zeros(len(fp), np.int8)
		#convert the fingerprint to a numpy array and wrap it into the dummy container
		DataStructs.ConvertToNumpyArray(fp, res)
		return FP(res)
	except:
		print("FPs for a structure cannot be calculated")
		return None

#Getting model/dataset stats for challenged set
def dists_yield(fps, nfps):
	# generator
	for i in range(1, nfps):
		yield [1-x for x in DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])]


#def ClusterData(fps, nPts, distThresh, isDistData=False, reordering=False):
def ClusterData(fps, nPts, distThresh, reordering=False):
	"""	clusters the data points passed in and returns the list of clusters

		**Arguments**

			- data: a list of items with the input data
				(see discussion of _isDistData_ argument for the exception)

			- nPts: the number of points to be used

			- distThresh: elements within this range of each other are considered
				to be neighbors
			- reodering: if this toggle is set, the number of neighbors is updated
					 for the unassigned molecules after a new cluster is created such
					 that always the molecule with the largest number of unassigned
					 neighbors is selected as the next cluster center.
		**Returns**
			- a tuple of tuples containing information about the clusters:
				 ( (cluster1_elem1, cluster1_elem2, ...),
					 (cluster2_elem1, cluster2_elem2, ...),
					 ...
				 )
				 The first element for each cluster is its centroid.

	"""
	nbrLists = [None] * nPts
	for i in range(nPts):
		nbrLists[i] = []

	#dmIdx = 0
	dist_fun = dists_yield(fps, nPts)
	for i in range(1, nPts):
		#print(i)
		dists = next(dist_fun)

		for j in range(i):
			#if not isDistData:
			#	dij = EuclideanDist(data[i], data[j])
			#else:
				#dij = data[dmIdx]
			dij = dists[j]
				#dmIdx += 1
			if dij <= distThresh:
				nbrLists[i].append(j)
				nbrLists[j].append(i)

	# sort by the number of neighbors:
	tLists = [(len(y), x) for x, y in enumerate(nbrLists)]
	tLists.sort(reverse=True)

	res = []
	seen = [0] * nPts
	while tLists:
		_, idx = tLists.pop(0)
		if seen[idx]:
			continue
		tRes = [idx]
		for nbr in nbrLists[idx]:
			if not seen[nbr]:
				tRes.append(nbr)
				seen[nbr] = 1
		# update the number of neighbors:
		# remove all members of the new cluster from the list of
		# neighbors and reorder the tLists
		res.append(tRes)
	return res



def ClusterFps(fps, method="Auto"):
	#Cluster size is probably smaller if the cut-off is larger. Changing its values between 0.4 and 0.45 makes a lot of difference
	nfps = len(fps)
	
	if method == "Auto":
		if nfps >= 10000:
			method = "TB"
		else:
			method = "Hierarchy"
	
	if method == "TB":
		#from rdkit.ML.Cluster import Butina
		cutoff = 0.56
		print("Butina clustering is selected. Dataset size is:", nfps)

		cs = ClusterData(fps, nfps, cutoff)
		
	elif method == "Hierarchy":
		print("Hierarchical clustering is selected. Dataset size is:", nfps)

		disArray = pdist(fps, 'jaccard')
		#Build model
		Z = hierarchy.linkage(disArray)
		
		#Cut-Tree to get clusters
		#x = hierarchy.cut_tree(Z,height = cutoff)
		average_cluster_size = 8
		cluster_amount = int( nfps / average_cluster_size )	 # calculate total amount of clusters by (number of compounds / average cluster size )
		x = hierarchy.cut_tree(Z, n_clusters = cluster_amount )		#use cluster amount as the parameter of this clustering algorithm. 
		
		#change the output format to mimic the output of Butina
		x = [e[0] for e in x]
		cs = [[] for _ in set(x)]

		for i in range(len(x)):
			cs[x[i]].append(i)
	return cs

#function for calculating similarities distributions between two datatsets
def SimHistogram(testSet, trainingSet):

	Max = []
	#for each compound in the test set, I will calculate the similarity to the nearest (most similar) training set compound
	for i in range(0, len(testSet)):
		sims = DataStructs.BulkTanimotoSimilarity(testSet[i], trainingSet)
		dists = sims
		dists.sort(reverse=True)
		Max.append(dists[0])

	#Build a distribution
	hs = np.histogram(Max, bins=10, range=(0.0, 1.0), density=False)

	return hs

def ChallengedSplit(data, smiles, fraction2train, clusterMethod="Auto", dropFPs=True, dropMol=True):
	data = data.copy()
	molecule = 'molecule'
	try:
		PandasTools.AddMoleculeColumnToFrame(data, smiles, molecule)
	except:
		print("Erroneous exception was raised and captured...")

	#remove records with empty molecules
	data = data.loc[data[molecule].notnull()]

	data['FP'] = [computeFP(m) for m in data[molecule]]

	#filter potentially failed fingerprint computations
	#data = data.loc[data['FP'].notnull()]

	fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in data[molecule]]
	if dropMol == True:
		data.drop(molecule, axis=1, inplace=True)	
	#data['FP'] = [computeFP(fp) for fp in fps]
	min_select = int(fraction2train * len(data))
	
	cluster_list = ClusterFps(fps, clusterMethod)

	cluster_list.sort(key=len, reverse=True)
	flat_list = sum(cluster_list, [])
	keep_tuple = flat_list[0 : min_select]
	
	if dropFPs == True:
		data.drop('FP', axis=1, inplace=True)

	train = data.iloc[list(keep_tuple)].copy()
	test = data.drop(train.index)
	Split = namedtuple('Split', 'testSet trainSet')
	return Split(testSet = test, trainSet = train)
	

#scale each dataframe row independently to pIC50 depending on unit type
class QueryAssay:
	"""
	# query by assay ID
	"""
	def __init__(self, mtable_dir):
		# initialize indices and MasterTable files
		self.index_file = '{}/PLS_Predictions/ZscaledAllAssaysAID_indices'.format(mtable_dir)
		self.MasterTable = '{}/PLS_Predictions/ZscaledAllAssaysAID.csv'.format(mtable_dir)
		
		self.indices = joblib.load(self.index_file)
		self.indices.index = self.indices.index.map(str)
		self.fid = open(self.MasterTable)
		self.separator = ',' if self.MasterTable.endswith('.csv') else '\t'

	def columns(self):
		self.fid.seek(0, 0)
		line = self.fid.readline().strip().split(self.separator)
		col = line[1: ]
		
		return(col)
		
	def idx(self):
		return(list(self.indices.index)[1: ])	

	def get(self, assay, raw=False):
		index = self.indices.loc[assay].values[0]
		self.fid.seek(index, 0)
		
		line = self.fid.readline().strip().split(self.separator)
		line_name = line[0]
		
		if raw:
			return(line_name, line[1: ])
		
		line_data = [float(x) for x in line[1: ]]
		
		return(line_name, line_data)


	
class QueryCmpd:
	"""
	# query by compound ID
	"""
	
	def __init__(self, mtable_dir):
		# initialize indices and MasterTable files
		self.index_file = '{}/{}'.format(mtable_dir, [fname for fname in os.listdir(mtable_dir) if fname.endswith('_indices')][0])
		self.MasterTable = '{}/{}'.format(mtable_dir, [fname for fname in os.listdir(mtable_dir) if fname.endswith('Table.csv')][0])
		
		self.indices = joblib.load(self.index_file)
		self.fid = open(self.MasterTable)
		self.separator = ',' if self.MasterTable.endswith('.csv') else '\t'

	def columns(self):
		self.fid.seek(0, 0)
		line = self.fid.readline().strip().split(self.separator)
		col = line[1: ]
		
		return(col)
		
	def get(self, cmpd, raw=False):
		index = self.indices.loc[cmpd].values[0]
		self.fid.seek(index, 0)
		
		line = self.fid.readline().strip().split(self.separator)
		line_name = line[0]
		
		if raw:
			return(line_name, line[1: ])
		
		line_data = [float(x) for x in line[1: ]]
		
		return(line_name, line_data)

class pls_dict:
	def __init__(self, model_list):
		self.dicts = self.make_dict(model_list)
	def make_dict(self, model_list):
		dicts = dict()
		for mdl in model_list:
			dicts[mdl] = joblib.load(mdl)
		return(dicts)
	def get_pls(self, mdl):
		model = self.dicts.get(mdl)
		return(model)
