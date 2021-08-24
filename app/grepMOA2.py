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


#!/usr/bin/env python
__doc__ = \
"""
assay/compound querying
========================================================================
run the following command on the zscaled file
Example (PRED_DIR is the location of ZscalledAllAssays)
#time grep -aob ".$" $PRED_DIR/ZscaledAllAssays.csv > $PRED_DIR/indices

"""

class QueryCmpd:
	"""
	input: mtable_dir, e.g., PLS_Predictions
	index file: MasterTable_indices
	RFR predicted master table: MasterTable.csv
	
	"""
	
	def __init__(self, mtable_dir, uncertainty=False):
		# initialize indices and MasterTable files
		if uncertainty == True:
			self.index_file = '{}/RFR_Predictions/ZscaledAllAssaysCID_indices'.format(mtable_dir)
			self.MasterTable = '{}/RFR_Predictions/ZscaledAllAssaysCID.csv'.format(mtable_dir)
		else:
			self.index_file = '{}/PLS_Predictions/ZscaledAllAssaysCID_indices'.format(mtable_dir)
			self.MasterTable = '{}/PLS_Predictions/ZscaledAllAssaysCID.csv'.format(mtable_dir)
		
		self.indices = joblib.load(self.index_file)
		self.fid = open(self.MasterTable)
		self.separator = ',' if self.MasterTable.endswith('.csv') else '\t'

	def columns(self):
		self.fid.seek(0, 0)
		line = self.fid.readline().strip().split(self.separator)
		col = line[1: ]
		return(col)
		
	def idx(self):
		return(list(self.indices.index)[1: ])
		
	def get(self, cmpd, raw=False):
		index = self.indices.loc[cmpd].values[0]
		self.fid.seek(index, 0)
		line = self.fid.readline().strip().split(self.separator)
		line_name = line[0]
		
		if raw:
			return(line_name, line[1: ])
		line_data = [float(x) if x != '' else 0.0 for x in line[1: ]]
		
		return(line_name, line_data)


class QueryAssay:
	"""
	input: mtable_dir, e.g., PLS_Predictions
	index file:	MasterTable_indices
	RFR predicted master table: MasterTable.csv
	
	"""
	
	def __init__(self, mtable_dir, uncertainty=False):
		# initialize indices and MasterTable files
		if uncertainty == True:
			self.index_file = '{}/RFR_Predictions/ZscaledAllAssaysAID_indices'.format(mtable_dir)
			self.MasterTable = '{}/RFR_Predictions/ZscaledAllAssaysAID.csv'.format(mtable_dir)
		else:
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

	

class QueryCustomCSV():
	def __init__(self, mtable, scale_stat, bool_ZpIC50):
		df_csv = pd.read_csv(mtable, header=0, index_col=0, sep=',')
		df_csv.index = df_csv.index.astype(str)
		df_stat = pd.read_csv(scale_stat, index_col=0, header=0, sep=',')

		if bool_ZpIC50 == True:
			cols = list(df_csv.columns)
			idx = [str(i) for i in df_stat.index]
			if len(set(cols) - set(idx)) == 0:
				cols = [int(c) for c in cols]
				df_mean = df_stat.loc[cols, 'mean_pred'].to_numpy()
				df_std = df_stat.loc[cols, 'stdev_pred'].to_numpy()
				df_csv = df_csv.sub(df_mean, axis=1).div(df_std, axis=1).round(3)
		self.df = df_csv

	def columns(self):
		col = list(self.df.columns)
		return(col)
	
	def get_column(self, assay=False, idx=False):
		if idx:
			return(list(self.df.index))
		else:
			return(self.df[assay])
		
	def get(self, cmpd):
		line_data = list(self.df.loc[cmpd])
		return(cmpd, line_data)
# test the compound's loading time

def get_CA(p, p_col, assays, cmpds):
	df_CA = pd.DataFrame(0.0, index=assays, columns=cmpds)
	cmpd_new = []
	
	for cmpd in cmpds:
		try:
			name, pqsar_vec = p.get(cmpd)
			df_get = pd.DataFrame([float(s) for s in pqsar_vec], index=p_col, columns=[name])
			df_CA[cmpd] = df_get.loc[assays]
		except:
			if cmpd in list(df_CA.columns):
				cmpd_new.append(cmpd)
			print('Warning! {} not found'.format(cmpd))
	
	df_CA.drop(cmpd_new, axis=1, inplace=True)

	return(df_CA)
		

def get_list(input_list):
	# get the query list
	separator = ',' if input_list.endswith('.csv') else '\t'
	df = pd.read_csv(input_list, header=0, index_col=0, sep=separator)
	items = [str(s) for s in set(df.index)]
	return(items)
	
def get_stat(df, thr, suffix='Zscore'):
	stats = pd.DataFrame(0.0, index=df.index, columns=[])
	col_name = list(df.columns)
	stats['count_{}>{}'.format(suffix, thr)] = df[df[col_name] >= thr].count(axis=1)
	stats['min_{}'.format(suffix)] = df.min(axis=1)
	stats['mean_{}'.format(suffix)] = df.mean(axis=1).round(3)
	stats['max_{}'.format(suffix)] = df.max(axis=1)
	
	return(stats)

def check_AID(items, p_col):
	real_items = [i for i in items if i in p_col]

	if len(real_items) == 0:
		print('Error! No AID was found')
		sys.exit(1)

	if len(real_items) < len(items):
		fake_items = [i for i in items if i not in real_items]
		print('Warning! AID {} not found'.format(fake_items))
	
	return(real_items)


def scale2pIC50(scale_stat, df, row_AID=True, uncertainty=False):
	df_scale = pd.read_csv(scale_stat, dtype={'AID':str}, index_col=0, header=0, sep=',')
	#df_scale.reindex([str(d) for d in df_scale.index])
	df_scale.index = df_scale.index.map(str)
	
	if row_AID:
		df = df.T.copy()

	cols = list(df.columns)
	df_std = df_scale.loc[cols, 'stdev_pred']
	if uncertainty == True:
		df = df * df_std
	else:
		df_mean = df_scale.loc[cols, 'mean_pred']
		df = df * df_std + df_mean

	if row_AID:
		df = df.T
	
	# output: cpds in the row, AID in the column
	return(df)

def tidy_view(df, args):
	# hide 3 columns: R^2_RF(ext), stdev_pred and mean_pred. I sort it vertically (descending) by Count_score>threshold and mean_score. 
	# I sort it horizontally (descending) by the column count_score>Threshold. 
	# If the columns are compounds, I calculate it in excel with =countif().

	cols = list(df.columns)
	for c in cols:
		if 'count_Zscore' in c or 'count_pIC50' in c:
			count_Zscore = c
		if 'mean_Zscore' in c or 'mean_pIC50' in c:
			mean_Zscore = c
	df.sort_values(by=[count_Zscore, mean_Zscore], axis=0, inplace=True, ascending=False)
	cpd_idx = []
	if args.Assay or args.CA:
		df.sort_values(by=[count_Zscore, mean_Zscore], axis=1, inplace=True, ascending=False)
		
		for idx in df.index:
			if 'count_Zscore' in idx or 'count_pIC50' in idx:
				break
			cpd_idx.append(idx)
	else:
		# args.Compound
		cpds = []
		for c in cols:
			if 'count_Zscore' in c or 'count_pIC50' in c:
				break
			cpds.append(c)
		AID = list(df.index)
		thr = float(count_Zscore.split('>')[1])

		df.loc['count_Zscore'] = 0.0
		df.loc['count_Zscore'] = df.loc[AID, cpds][df.loc[AID, cpds] >= thr].count()
		df.loc['count_Zscore'] = [float(s) for s in df.loc['count_Zscore']]
		df = df.sort_values(by='count_Zscore', axis=1, ascending=False)
		df.drop('count_Zscore', axis=0, inplace=True)
		#tobe_moved = ['R^2_RF(ext)', 'stdev_pred', 'mean_pred', 'validity', 'delta_error_rate', 'efficiency', 'wt']
		tobe_moved = ['stdev_pred', 'mean_pred']
		cols = list(df.columns)
		left_cols = [c for c in cols if c not in tobe_moved]
		left_cols.extend(tobe_moved)
		df = df[left_cols]
	
	return(df, cpd_idx)


def save_query(df, args, out_csv, scale_stat, bool_ZpIC50):

	cols = []
	rows = []
	for c in df.columns:
		if 'count_' in c:
			break
		cols.append(c)
	
	cols_left = [c for c in df.columns if c not in cols]

	for i in df.index:
		if 'count_' in i:
			break
		rows.append(i)
	
	if args.Uncertainty.lower() in ['true', 'yes', 't', 'y']:
		if args.Local:
			local_dir, basename = os.path.dirname(args.Local), os.path.basename(args.Local)
			local_csv_PI = os.path.join(local_dir, basename.split(',')[1])
			u = QueryCustomCSV(local_csv_PI, scale_stat, False)
			u_col = u.columns()
		else:
			u = QueryCmpd(args.Directory, uncertainty=True)
			u_col = u.columns()


		df_error = get_CA(u, u_col, rows, cols) if args.Compound else get_CA(u, u_col, cols, rows)
		df_error = scale2pIC50(scale_stat, df_error, uncertainty=True) if bool_ZpIC50 == False else df_error

		df_error = df_error if args.Compound else df_error.T

		error_cols = list(df_error.columns)
		error_cols = [col + '_Error' for col in error_cols]
		df_error.columns = error_cols

		df = df.merge(df_error, how='left', right_index=True, left_index=True)
		
		assert len(cols) == len(error_cols)
		cols_all = sum([[c, e] for c, e in zip(cols, error_cols)], []) + cols_left
		df = df[cols_all]

	if args.Compound:
		df.round(2).to_csv(out_csv, encoding='utf-8')
	else:
		df.round(2).loc[rows].to_csv(out_csv, mode='a+', encoding='utf-8')
		with open(out_csv, 'a+') as fid:
			fid.writelines('\n')
		df.drop(rows, axis=0).round(2).to_csv(out_csv, mode='a+', header=False, encoding='utf-8')

def main():
	description = """Querying individual screening of pQSAR by CID and/or AID
	***usage examples
	***input file (.txt or .csv): querying by compounds or assays, a header line followed by compounds or assays (one item per line)
				querying by compounds and assays, a header line followed by CID at the first column and AID second column. 

	***querying by compounds (CID) with a threshold of 3 stdev above the mean
	python grepMOA2.py -c -i cid.txt -d ../chembl_28 -t 3.0 -z True -o cid_out.csv
	
	***querying by assays (AID) with a threshold of 3 stdev above the mean
	python grepMOA2.py -a -i aid.txt -d ../chembl_28 -t 3.0 -z True -o aid_out.csv

	***querying by compounds (CID) and assays (AID)
	python grepMOA2.py -ca -i ca_id.txt -d ../chembl_28 -z True -o ca_id_out.csv
"""
	
	epilog = """----------------------profile-QSAR application-----------------------
		"""
	parser = argparse.ArgumentParser(
						description=description,
						formatter_class=argparse.RawDescriptionHelpFormatter,
						epilog=epilog)


	parser.add_argument('-i', '--Input', action='store', help='Input file with header line followed by querying compounds\' ID and/or assays\' AIDs', metavar='')
	parser.add_argument('-d', '--Directory', help='Directory contains modeling information, default is chembl_28', type=str, metavar='')
	parser.add_argument('-u', '--Uncertainty', help='Query the predicted uncertainty (Error) at 80 percent confidence level', type=str, default='False', metavar='True/False')
	parser.add_argument('-l', '--Local', help='Local file (csv) of custom prediction', type=str, metavar='')
	parser.add_argument('-e', '--Experimental', help='False(default): Querying predicted values; True: Querying experimental values (always returns pIC50, -z is not applicable)', type=str, default='False', metavar='True/False')
	
	#parser.add_argument('-z', '--ZpIC50', help='Z scaled predictions of real ones', action='store_true')

	parser.add_argument('-z', '--ZpIC50', help='True (default): Threshold and predictions in Z-scaled values; False: original pIC50 (log molar).', type=str, default='True', metavar='True/False')
	parser.add_argument('-t', '--Threshold', help='Threshold to filter out unqualified screening', metavar='')
	parser.add_argument('-o', '--Output', help='Output file in csv format', default='Query_Output.csv', metavar='')
	
	index_group = parser.add_mutually_exclusive_group()
	index_group.add_argument('-c', '--Compound', action='store_true', help='Querying by compounds (one compound per row)')
	index_group.add_argument('-a', '--Assay', action='store_true', help='Querying by assays (one assay per row)')
	index_group.add_argument('-ca', '--CA', action='store_true', help='Querying by compounds(first column) and assays (second column)')

	if len(sys.argv) < 4:
		parser.print_help()
		sys.exit(1)

	args = parser.parse_args()
	Uncertainty = True if args.Uncertainty.lower() in ['true', 'yes', 't', 'y'] else False
	if not args.Input:
		print('Error: no input file ')
		sys.exit(1)
	

	# Default threshold is 3.0
	#Thr = args.Threshold if args.Threshold else 3.0
	# flag predicted uncertainty (error)
	# Default output
	
	bool_expr = True if args.Experimental.lower() in ['true', 'yes', 't', 'y'] else False
	out_csv = args.Output if args.Output else 'Query_Output.csv'
	if os.path.exists(out_csv):
		os.remove(out_csv)
	
	if bool_expr == False:
	
		bool_ZpIC50 = True if args.ZpIC50.lower() in ['true', 'yes', 't', 'y'] else False
		scale_stat = '{}/PLS_Predictions/ZscalingStats.csv'.format(args.Directory)
		summary_table = '{}/Summary_pQSAR/SummaryPQSAR.csv'.format(args.Directory)
		df_summary = pd.read_csv(summary_table, index_col=0, header=0)
		df_summary.index = df_summary.index.map(str)

		if args.Local:
			if Uncertainty == True:
				# we expect AllPLS_prediction.csv and AllPLS_prediction_PI80.csv are comma seperated.
				# e.g., AllPLS_prediction.csv,AllPLS_prediction_PI80.csv
				local_dir, basename = os.path.dirname(args.Local), os.path.basename(args.Local)
				local_csv_act = os.path.join(local_dir, basename.split(',')[0])
			else:
				local_csv_act = args.Local

			if local_csv_act.endswith('.csv'):
				p = QueryCustomCSV(local_csv_act, scale_stat, bool_ZpIC50)
				p_col = p.columns()
			else:
				print('Error! Need to be a custom predicted csv file')
				sys.exit(1)
			if args.Assay:
				p_idx = p.get_column(idx=True)
		else:
			if args.Compound or args.CA:
				p = QueryCmpd(args.Directory)
				p_col = p.columns()
			elif args.Assay: 
				p = QueryAssay(args.Directory)
				p_idx = p.columns()
				p_col = p.idx()

		# get the list of compounds or assays
		items = get_list(args.Input)

		if args.Compound:
			items = [it.strip() for it in items]
			df_cmpd = pd.DataFrame(0.0, index=p_col, columns=items)
			cmpd_new = []

			for cmpd in items:
				#print(cmpd, thr)
				try:
					name, pqsar_vec = p.get(cmpd)
					df_cmpd[cmpd] = pqsar_vec
				except:
					cmpd_new.append(cmpd)
					print('Warning! {} not found'.format(cmpd))
			df_cmpd.drop(cmpd_new, axis=1, inplace=True)
			
			if len(df_cmpd.index) == 0:
				print('Error, No CID was found')
				sys.exit(1)

			if not args.Local and bool_ZpIC50 == False:
					# convert back to pIC50
				df_cmpd = scale2pIC50(scale_stat, df_cmpd)

			if not args.Threshold:
				thr= -10000
			else:
				thr = float(args.Threshold)
				df_cmpd = df_cmpd.iloc[list(df_cmpd.max(axis=1) > thr), :].copy()

			# If threshold argument is not provided, the whole result table will be saved by default.
			if bool_ZpIC50 == False:
				stats = get_stat(df_cmpd, thr, suffix='pIC50')
			else:
				#print(df_cmpd.iloc[0:5, 0:10])
				stats = get_stat(df_cmpd, thr)

			df_cmpd = df_cmpd.merge(stats, how='left', left_index=True, right_index=True)
			
			df_cmpd.index = df_cmpd.index.map(str)
			df_cmpd = df_cmpd.merge(df_summary, how='left', left_index=True, right_index=True)
			df_cmpd.index.name = 'AID'
			df_cmpd, _ = tidy_view(df_cmpd, args)
			save_query(df_cmpd, args, out_csv, scale_stat, bool_ZpIC50)

		if args.Assay:
			
			items = check_AID(items, p_col)
			df_assay = pd.DataFrame(0.0, index=p_idx, columns=items)

			# filtering z scaled values
			if args.Local and os.path.isfile(local_csv_act):
				df_assay[items] = p.get_column([str(s) for s in items])
			
			else:
				for assay in items:
					name, pqsar_vec = p.get(str(assay))
					df_assay[assay] = pqsar_vec
		
				if bool_ZpIC50 == False and not args.Local:
				# convert back to pIC50
					df_assay = scale2pIC50(scale_stat, df_assay, row_AID=False)

			
			#df_zscore = (df_assay - df_assay.mean(axis=0)) / df_assay.std(axis=0)

			thr = float(args.Threshold) if args.Threshold else -1000	
			df_assay_thr = df_assay.iloc[list(df_assay.max(axis=1) > thr), :]
			cpd_idx = list(df_assay_thr.index)
			
			if bool_ZpIC50 == False:
				stats_assay = get_stat(df_assay_thr, thr, suffix='pIC50')
			else:
				stats_assay = get_stat(df_assay_thr, thr)
			
			df_assay_thr = df_assay_thr.T.copy()
			
			if bool_ZpIC50 == False:
				stats_cmpd = get_stat(df_assay_thr, thr, suffix='pIC50')
			else:
				stats_cmpd = get_stat(df_assay_thr, thr)
				
			df_assay_thr = df_assay_thr.merge(stats_cmpd, how='left', left_index=True, right_index=True)
			df_assay_thr = df_assay_thr.merge(df_summary, how='left', left_index=True, right_index=True)
			df_assay_thr = df_assay_thr.T.copy()
			df_assay_thr = df_assay_thr.merge(stats_assay, how='left', left_index=True, right_index=True)
			
			df_assay_thr.index.name = 'CID'
			df_assay_thr, cpd_idx = tidy_view(df_assay_thr, args)

			save_query(df_assay_thr, args, out_csv, scale_stat, bool_ZpIC50)
			

		if args.CA:
			# the first column is CID and the second column AIDs

			thr = float(args.Threshold) if args.Threshold else 3.0
			cmpds = []
			assays = []
			separator = ',' if args.Input.endswith('.csv') else '\t'
			with open(args.Input) as fid:
			#with open(input_list) as fid:
				fid.readline()
				
				for line in fid:
					if not line.startswith(separator):
						lf = line.strip().split(separator)
						
						if lf[0] != '':
							cmpds.append(lf[0].strip())
						if len(lf) > 1 and lf[1] != '':
							assays.append(lf[1])
					else:
						lf = line.strip().strip(separator)
						assays.append(lf)
			
			assays = check_AID(assays, p_col)
			df_get = get_CA(p, p_col, assays, cmpds)
			
			if not args.Local and bool_ZpIC50 == False:
				df_get = scale2pIC50(scale_stat, df_get)

			if bool_ZpIC50 == False:
				stats_assay = get_stat(df_get, thr, suffix='pIC50')
			else:
				stats_assay = get_stat(df_get, thr)

			cpd_idx = list(df_get.columns)
			df_get = df_get.T.copy()
			
			if bool_ZpIC50 == False:
				stats_cmpd = get_stat(df_get, thr, suffix='pIC50')
			else:
				stats_cmpd = get_stat(df_get, thr)	
			
			df_get = df_get.merge(stats_cmpd, how='left', left_index=True, right_index=True)
			df_get = df_get.T.copy()
			df_get = df_get.merge(stats_assay, how='left', left_index=True, right_index=True)
			df_get = df_get.merge(df_summary, how='left', left_index=True, right_index=True)
			df_get = df_get.T.copy()
			df_get.index.name = 'CID'
			df_get, cpd_idx = tidy_view(df_get, args)

			save_query(df_get, args, out_csv, scale_stat, bool_ZpIC50)

	else:
		###################################################
		"""
		Experimental pIC50 querying.
		"""		
		##################################################
		summary_table = '{}/Summary_pQSAR/SummaryPQSAR.csv'.format(args.Directory)
			#df_summary = pd.read_excel(summary_table, index_col=0, header=0)
		df_summary = pd.read_csv(summary_table, index_col=0, header=0)
		df_summary.index = df_summary.index.map(str)
		
		## only work for version chembl_28 in this particular case
		file_path = '{}/chembl_28_AllAssayData.txt'.format(args.Directory)


		df_expr = pd.read_csv(file_path, index_col=None, header=0, sep='\t')
		col_names = list(df_expr.columns)
		col_names[0: 3] = ['CID', 'AID', 'PIC50']
		df_expr['AID'] = df_expr['AID'].astype(str)
		df_expr.columns = col_names
		df_expr = df_expr[col_names[0 : 3]]
		thr = float(args.Threshold) if args.Threshold else -10000
		
		items = get_list(args.Input)
		if args.Compound:
			df_cmpd = df_expr[df_expr['CID'].isin(items)]
			
			if args.Threshold:
				df_cmpd_thr = df_cmpd[df_cmpd['PIC50'] > thr]
			else:
				df_cmpd_thr = df_cmpd
			
			df_cmpd_thr = pd.pivot_table(df_cmpd_thr, index=df_cmpd_thr['AID'], columns=df_cmpd_thr['CID'])
			df_cmpd_thr.columns = list(df_cmpd_thr.columns.levels[1])
			
			stats = get_stat(df_cmpd_thr, -10000, suffix='pIC50')
			df_cmpd_thr = df_cmpd_thr.merge(stats, how='left', left_index=True, right_index=True)
			df_cmpd_thr = df_cmpd_thr.merge(df_summary, how='left', left_index=True, right_index=True)
			df_cmpd_thr.index.name = 'AID'

			df_cmpd_thr, _ = tidy_view(df_cmpd_thr, args)
			df_cmpd_thr.to_csv(out_csv, encoding='utf-8')		

		if args.Assay:
			df_assay = df_expr[df_expr['AID'].isin(items)]
			
			if args.Threshold:
				df_assay_thr = df_assay[df_assay['PIC50'] > thr]
			else:
				df_assay_thr = df_assay
			
			df_assay_thr = pd.pivot_table(df_assay_thr, index=df_assay_thr['CID'], columns=df_assay_thr['AID'])
			df_assay_thr.columns = list(df_assay_thr.columns.levels[1])
			cpd_idx = list(df_assay_thr.index)
			stats_assay = get_stat(df_assay_thr, thr, suffix='pIC50')
			df_assay_thr = df_assay_thr.T.copy()

			stats_cmpd = get_stat(df_assay_thr, thr, suffix='pIC50')
			
			df_assay_thr = df_assay_thr.merge(stats_cmpd, how='left', left_index=True, right_index=True)
			df_assay_thr = df_assay_thr.merge(df_summary, how='left', left_index=True, right_index=True)
			df_assay_thr = df_assay_thr.T.copy()
			df_assay_thr = df_assay_thr.merge(stats_assay, how='left', left_index=True, right_index=True)
			
			
			df_assay_thr.index.name = 'CID'
			df_assay_thr, cpd_idx = tidy_view(df_assay_thr, args)
			df_assay_thr.loc[cpd_idx].to_csv(out_csv, mode='a+', encoding='utf-8')
			with open(out_csv, 'a+') as fid:
				fid.writelines('\n')
			df_assay_thr.drop(cpd_idx, axis=0).to_csv(out_csv, mode='a+', header=False, encoding='utf-8')
			
		if args.CA:
			separator = ',' if args.Input.endswith('.csv') else '\t'
			cmpds = []
			assays = []
			
			with open(args.Input) as fid:
			#with open(input_list) as fid:
				fid.readline()
				
				for line in fid:
					if not line.startswith(separator):
						lf = line.strip().split(separator)
						
						if lf[0] != '':
							cmpds.append(lf[0])
						if len(lf) > 1 and lf[1] != '':
							assays.append(lf[1])
					else:
						lf = line.strip().strip(separator)
						assays.append(lf)
		
			df_expr = df_expr[df_expr['AID'].isin(assays)]
			df_expr = df_expr[df_expr['CID'].isin(cmpds)]
			
			
			if args.Threshold:
				df_expr_thr = df_expr[df_expr['PIC50'] > thr]
			else:
				df_expr_thr = df_expr
			
			df_expr_thr = pd.pivot_table(df_expr_thr, index=df_expr_thr['CID'], columns=df_expr_thr['AID'])
			df_expr_thr.columns = list(df_expr_thr.columns.levels[1])
			cpd_idx = list(df_expr_thr.index)
			stats_assay = get_stat(df_expr_thr, thr, suffix='pIC50')
			df_expr_thr = df_expr_thr.T.copy()

			stats_cmpd = get_stat(df_expr_thr, thr, suffix='pIC50')
			
			df_expr_thr = df_expr_thr.merge(stats_cmpd, how='left', left_index=True, right_index=True)
			df_expr_thr = df_expr_thr.merge(df_summary, how='left', left_index=True, right_index=True)
			df_expr_thr = df_expr_thr.T.copy()
			df_expr_thr = df_expr_thr.merge(stats_assay, how='left', left_index=True, right_index=True)
			
			df_expr_thr.index.name = 'CID'

			df_expr_thr, cpd_idx = tidy_view(df_expr_thr, args)

			df_expr_thr.loc[cpd_idx].to_csv(out_csv, mode='a+', encoding='utf-8')
			with open(out_csv, 'a+') as fid:
				fid.writelines('\n')
			df_expr_thr.drop(cpd_idx, axis=0).to_csv(out_csv, mode='a+', header=False, encoding='utf-8')
		
if __name__ == "__main__":
	import os, sys
	import time
	import argparse
	import subprocess
	import pandas as pd
	import joblib
	main()
