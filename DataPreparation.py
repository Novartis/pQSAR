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
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
#from rdkit.Chem.MolStandardize import rdMolStandardize

def fragment_remover(smi):
	"""
	# deal with molecule with more than one covalently bonded unit 
	# smiles separated by dot (*.***):
	# 		NCCO.NC(=O)c1ccc2cc(CCC3=C(NC(=O)NC3=O)C(=O)O)ccc2c1
	# 		[Na+].COc1nc(cc2nc([n-]c12)c3c(Cl)c(nn3C)C(C)(C)C)c4c(Cl)cccc4Cl
	# remove salts etc. 
	# keep the largest covalent unit in a molecule with multiple fragments
	# 
	"""
	# only compare the length of smiles of fragments
	# randomly pick one if the length of two fragments(smiles) equals

	if smi.find('.') != -1:
		smi_frag = smi.split('.')
		#smi_frag = [it.strip('+]') for it in smi_frag]
		smi = max(smi_frag, key = len)
	
	smi = smi.replace('@', '')
	return(smi)

def smiles_validity(cmpd):
	cmpd_len = len(cmpd)
	df = []
	flag = True
	count = 0
	chunk = 50_000

	while flag:
		if count + chunk < cmpd_len:
			tmp = cmpd.iloc[range(count, count + chunk)].copy()
			count += chunk
		else:
			tmp = cmpd.iloc[count:,].copy()
			flag = False
		try:
			PandasTools.AddMoleculeColumnToFrame(tmp, 'smi', 'molecule')
		except:
			print("Erroneous exception was raised and captured...")
		
		tmp = tmp.loc[tmp['molecule'].notnull()]
		tmp['SMILES'] = tmp.apply(lambda x: Chem.MolToSmiles(x[3]), axis=1)
		tmp.drop(['smi', 'molecule'], axis=1, inplace=True)
		df.append(tmp)
	df = pd.concat(df)
	
	return df

cmpd_file = sys.argv[1]
act_file = sys.argv[2]
cmpd_dir = sys.argv[3]
act_dir = sys.argv[4]
if not os.path.isfile(cmpd_dir):
	cmpd = pd.read_csv(cmpd_file, header=0, index_col=None)
	col_names = list(cmpd.columns)
	col_names[0: 2] = ['PREFERRED_NAME', 'SMILES']
	cmpd.columns = col_names
	cmpd = cmpd[col_names[0: 2]]
	cmpd.dropna(axis=0, inplace=True)
	cmpd['smi'] = cmpd.apply(lambda x: fragment_remover(x[1]), axis=1)
	cmpd = smiles_validity(cmpd)

	cmpd.to_csv(cmpd_dir, index=False, header=True)
else:
	cmpd = pd.read_csv(cmpd_dir, header=0, index_col=None)

# CID  AID  PIC50   smiles  Units   ResultType  QUALIFIER
if not os.path.isfile(act_dir):
	sept = ',' if act_file.endswith('.csv') else '\t'
	data = pd.read_csv(act_file, index_col=None, header=0, sep=sept)
	col_name = list(data.columns)
	col_name[0: 4] = ['CID', 'AID', 'PIC50', 'smiles']
	data.columns = col_name
	data = data.merge(cmpd, left_on='CID', right_on='PREFERRED_NAME')
	data['smiles'] = data['SMILES']
	data.drop(['SMILES', 'PREFERRED_NAME'], axis=1, inplace=True)
	data.dropna(axis=0, how='any', subset=['CID', 'AID', 'PIC50', 'smiles'], inplace=True)
	data = data.sort_values(by=['AID', 'smiles', 'PIC50'], axis=0).drop_duplicates(subset=['AID', 'smiles'], keep='last')
	data.round(3).to_csv(act_dir, index=False, header=True, sep=sept)
