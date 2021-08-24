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


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import sys

import numpy as np
import pandas as pd
import joblib
import CommonTools as Tools

assay_file = sys.argv[1]
mdl_addr = sys.argv[2]
mdl_type = sys.argv[3]
file2predict = sys.argv[4]
master_addr = sys.argv[5]
index_save = bool(sys.argv[6])
NVSReference = sys.argv[7]

if mdl_type == 'RF':
	master_table = '{}/{}_master_table_tmp.csv'.format(master_addr, assay_file.split('/')[-1])
	
	if os.path.exists(master_table):
		sys.exit()
	
	assay_list = [s.strip().split('\t')[0] for s in open(assay_file)]
	model_list = ['{}/{}.RFmodel'.format(mdl_addr, a) for a in assay_list]
	df_NVS = pd.read_csv(NVSReference, index_col=0)
	df_master_table = pd.DataFrame(0.0, index=df_NVS.index, columns=assay_list)
	rf_models = Tools.pls_dict(model_list)
	
	if ',' in file2predict:
		file2predict = file2predict.split(',')
	else:
		file2predict = [file2predict]
	
	for fps in file2predict:
		cmpd_fps = joblib.load(fps)
		cmpd_idx = cmpd_fps.index

		for k, mdl in enumerate(model_list):
			model = rf_models.get_pls(mdl)

			data_tmp = pd.read_csv('{}/{}_tmp_split.txt'.format(mdl_addr, assay_list[k]), header=0, index_col=0, sep='\t')
			data_tmp = data_tmp[~data_tmp.index.duplicated()]
			col_names = list(data_tmp.columns)
			col_names[0 : 3] = ['AssayID', 'PIC50', 'smiles']
			idx_tmp = set(data_tmp.index)

			data_tmp.columns = col_names
			df_master_table.loc[cmpd_idx, assay_list[k]] = model.predict(cmpd_fps).astype('float32')
			comm_idx = idx_tmp.intersection(set(df_master_table.index))
			if len(comm_idx) > 0:
				# replace RF fit with experimental measurements
				df_master_table.loc[comm_idx, assay_list[k]] = data_tmp.loc[comm_idx, 'PIC50']

	del cmpd_fps
	del model
	"""
	master table storage
	"""
	#df_master_table.drop('CID', axis=1, inplace=True)
	
	df_master_table.sort_index(inplace=True)
	
	df_master_table.to_csv(master_table, index=False, float_format='%.3f')

elif mdl_type == 'PLS':
	chunkSize = 50_000
	##TEST:
	#chunkSize = 1000
	master_tableNVS = '{}/{}_master_table_tmp_CID.csv'.format(master_addr, assay_file.split('/')[-1])
	master_tableCDS = '{}/{}_master_table_tmp_AID.csv'.format(master_addr, assay_file.split('/')[-1])
	stat_table = '{}/{}_stat_table_tmp.csv'.format(master_addr, assay_file.split('/')[-1])
	
	if os.path.exists(master_tableNVS) and os.path.exists(master_tableCDS):
		sys.exit()

	# file2predict here should be the master table
	assay_list = [s.strip().split('\t')[0] for s in open(assay_file)]
	

	model_list = ['{}/{}.PLSmodel'.format(mdl_addr, a) for a in assay_list]

	pls_models = Tools.pls_dict(model_list)
	
	df_pilot_3 = pd.read_csv(file2predict, header=0, index_col=0, nrows=1)
	float_cols = [c for c in df_pilot_3 if df_pilot_3[c].dtype == "float64"]
	float32_cols = {c: np.float32 for c in float_cols}
	
	del df_pilot_3
	
	iterate = pd.read_csv(file2predict, header=0, index_col=0, engine='c', dtype=float32_cols, iterator=True, chunksize=chunkSize)

	df_master_table = []
	
	for i, mat in enumerate(iterate):
		tmp = pd.DataFrame(0.0, index=mat.index, columns=assay_list, dtype=np.float32)
		
		for k, mdl in enumerate(model_list):
			model = pls_models.get_pls(mdl)
			assay = assay_list[k]
			dscr = model[-2]
			if 'compoundID' in dscr:
				dscr.remove('compoundID')
			
			if 'PIC50' in dscr:
				dscr.remove('PIC50')
			if 'smiles' in dscr:
				dscr.remove('smiles')		
			
			model = model[0]
			dscr = [d.strip() for d in dscr]
			mat_dscr = mat[dscr].copy()
			tmp[assay] = [m[0] for m in model.predict(mat_dscr).astype('float32')]

		df_master_table.append(tmp)

	del mat
	del tmp
	del mat_dscr
	del iterate
	del pls_models
	
	df_master_table = pd.concat(df_master_table, axis=0)

	STD = df_master_table.std()
	MEAN = df_master_table.mean()

	
	df_stat = pd.DataFrame(0.0, index=df_master_table.columns, columns=['stdev_pred', 'mean_pred', 'max_pred', 'min_pred', 'Q99.5', 'Q99', 'Q95', 'Q75', 'Q25', 'Q05', 'Q01', 'Q0.5'], dtype=np.float32)
	
	df_stat['stdev_pred'] = STD
	df_stat['mean_pred'] = MEAN
	df_stat['max_pred'] = df_master_table.max()
	df_stat['min_pred'] = df_master_table.min()
	
	df_stat['Q99.5'] = df_master_table.quantile(0.995)
	df_stat['Q99'] = df_master_table.quantile(0.99)
	df_stat['Q95'] = df_master_table.quantile(0.95)
	df_stat['Q75'] = df_master_table.quantile(0.75)
	df_stat['Q25'] = df_master_table.quantile(0.25)
	df_stat['Q05'] = df_master_table.quantile(0.05)	
	df_stat['Q01'] = df_master_table.quantile(0.01)
	df_stat['Q0.5'] = df_master_table.quantile(0.005)
	
	df_stat.index.name = 'AssayID'
	df_master_table = (df_master_table - MEAN) / STD
	
	df_master_table.sort_index(inplace=True)
	df_master_table.to_csv(master_tableNVS, index=False, float_format="%.3f")
	df_master_table = df_master_table.T
	df_master_table.index.name = 'AssayID'
	df_master_table.to_csv(master_tableCDS, index=True, float_format="%.3f")
	df_stat.to_csv(stat_table, index=True, float_format="%.3f")
