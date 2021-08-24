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
"""
#this is python 3.8.6!
"""
#######################################################
## classes and functions

def get_job_id(process):
	job_sub = process.stdout.decode('utf-8')
	try: 
		jobID = re.search('Your job-array (\d.+?)[ \.]', job_sub)[1]
	except:
		jobID = re.search('Your job (\d.+?)[ \.]', job_sub)[1]

	return(jobID)


def paste_pred(resultsDir, ResultExtension, CIDReference, MasterTable, Job_name):
	# paste all the predict file horizontally into one
	# first use 26 cluster jobs and then just one job
	merge_round_one = '_Cmpd_merge_round_one_tmp'
	merge_round_two = '_Cmpd_merge_round_two_tmp'
	prefixes = list(set([f[0 : 6] for f in os.listdir(resultsDir) if f.endswith(ResultExtension)]))
	jobIDs = []

	for prefix in prefixes:
		bash_line = get_bash(job_name=Job_name + '1', job_mem='8G')
		#bash_line += "awk " + "'" + "{ a[FNR] = (a[FNR] ? a[FNR] FS : \"\") $2 } END { for(i=1;i<=FNR;i++) print a[i] }" + "'" + " {}/{}*pred_tmp.txt > {}/merge_round_one_{}_tmp.txt\n".format(tmp_dir, i, tmp_dir, i)
		bash_line += 'cd {}\n'.format(resultsDir)
		bash_line += 'paste -d , {}*{} > {}{}.csv\n'.format(prefix, ResultExtension, prefix, merge_round_one)
		
		bash_file = '{}/{}{}.sh'.format(resultsDir, prefix, merge_round_one)
		with open(bash_file, 'w') as fid:
			fid.writelines(bash_line)

		process = subprocess.run(["qsub", "-l", "m_mem_free=" + '8G', bash_file], stdout=subprocess.PIPE)
		jobID = get_job_id(process)
		jobIDs.append(jobID)

	check_current_jobs([jobID])
	##TEST:
	#check_current_jobs([jobID], sec=5)

	bash_file = 'merge_round_two_tmp.sh'
	bash_line = get_bash(job_name=Job_name + '2', job_mem='8G')
	bash_line += 'cd {}\n'.format(resultsDir)
	bash_line += "paste -d , {} {}/*{}.csv > {}/{}\n".format(CIDReference, resultsDir, merge_round_one, resultsDir, MasterTable)
	bash_file = '{}/{}{}.sh'.format(resultsDir, Job_name, merge_round_two)
	with open(bash_file, 'w') as fid:
		fid.writelines(bash_line)

	process = subprocess.run(["qsub", "-l", "m_mem_free=8G", bash_file], stdout=subprocess.PIPE)
	
	jobID = get_job_id(process)
	check_current_jobs([jobID])
	##TEST:
	#check_current_jobs([jobID], sec=5)

	
def paste_pred2(resultsDir, ResultExtension, MasterTable, Job_name):
	# paste all the predict file vertically into one
	# first use 26 cluster jobs and then just one job
	merge_round_one = '_AID_PLS_merge_round_one_tmp'
	merge_round_two = '_AID_PLS_merge_round_two_tmp'
	
	prefixes = list(set([f[0 : 6] for f in os.listdir(resultsDir) if f.endswith(ResultExtension)]))
	#jobIDs = []

	for prefix in prefixes:
		bash_line = get_bash(job_name=Job_name + '1', job_mem='8G')
		bash_line += 'cd {}\n'.format(resultsDir)
		bash_line += 'RDN=`ls {}*{}| head -1`\n'.format(prefix, ResultExtension)
		bash_line += 'head -1 $RDN > {}{}\n'.format(prefix, merge_round_one)
		bash_line += 'C5=`head -1 $RDN | cut -c1-5`\n'
		bash_line += "cat {}*{} | grep -v $C5 >> {}{}\n".format(prefix, ResultExtension, prefix, merge_round_one)

		bash_file = '{}/{}{}.sh'.format(resultsDir, prefix, merge_round_one)
		with open(bash_file, 'w') as fid:
			fid.writelines(bash_line)

		process = subprocess.run(["qsub", "-l", "m_mem_free=" + '8G', bash_file], stdout=subprocess.PIPE)
		jobID = get_job_id(process)

	check_current_jobs([jobID])
	##TEST:
	#check_current_jobs([jobID], sec=5)

	bash_file = 'merge_round_two_tmp.sh'
	bash_line = get_bash(job_name=Job_name + '2', job_mem='8G')
	bash_line += 'cd {}\n'.format(resultsDir)
	bash_line += 'RDN=`ls *{}| head -1`\n'.format(merge_round_one)
	bash_line += 'head -1 $RDN > {}\n'.format(MasterTable)
	bash_line += 'C5=`head -1 $RDN | cut -c1-5`\n'
	bash_line += 'cat *{} | grep -v $C5 >> {}\n'.format(merge_round_one, MasterTable)	
	bash_file = '{}/{}{}.sh'.format(resultsDir, Job_name, merge_round_two)
	
	with open(bash_file, 'w') as fid:
		fid.writelines(bash_line)

	process = subprocess.run(["qsub", "-l", "m_mem_free=8G", bash_file], stdout=subprocess.PIPE)
	
	jobID = get_job_id(process)
	check_current_jobs([jobID])
	##TEST:
	#check_current_jobs([jobID], sec=5)

def check_current_jobs(jobIDs, sec=60):
	# jobIDs is a list
	while True:
		jobs = subprocess.Popen("qstat | tail -n +3 | awk -F' ' '{print $1}' | sort -u", shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')
		jobs.remove('')
		running = subprocess.Popen("qstat | tail -n +3 | awk -F' ' '{print $5}' | sort -u", shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')
		running.remove('')
		
		if len(set(jobIDs).intersection(set(jobs))) == 0:
			break
		
		if len(running) == 1 and running[0] == 'dr':
			break
		else:
			time.sleep(sec)

def get_bash(job_name, job_num=1, job_mem='4G', job_time=340000, threads=1, py=False):
	# generate the first few lines of bash file 
	bash_line = '#$ -S /bin/bash\n'
	bash_line += '#$ -N {}_{}\n'.format(job_name, job_num)
	bash_line += '#$ -cwd\n'
	bash_line += '#$ -o o_{}_tmp.txt\n'.format(job_name)	
	bash_line += '#$ -e e_{}_tmp.txt\n'.format(job_name)
	bash_line += '#$ -q default.q\n'
	bash_line += '#$ -l h_rt={}\n'.format(job_time)
	bash_line += '#$ -j y\n'
	bash_line += '#$ -V\n'
	bash_line += '#$ -t 1:{}\n'.format(job_num)
	#bash_line += '#$ -tc 200\n'
	bash_line += '#$ -l m_mem_free={}\n'.format(job_mem)
	bash_line += '#$ -pe smp {}\n'.format(threads)
	#bash_line += '#$ -binding linear:{}\n\n'.format(threads)
	bash_line += 'export OMP_NUM_THREADS={}\n'.format(threads)
	bash_line += 'export OPENBLAS_NUM_THREADS={}\n'.format(threads)
	bash_line += 'export MKL_NUM_THREADS={}\n\n'.format(threads)

	return(bash_line)

	
def split_by_cmpd_number(df_info, CMPD_LIMIT, MDL_mem, out_dir, cmpd_num, threads=1):
	
	assay_class = []
	assay_num_in_class = []
	
	for i, num in enumerate(CMPD_LIMIT):
		if i == 0:
			df_tmp = df_info.loc[df_info[cmpd_num] >= CMPD_LIMIT[i]]
			assay_class.append('assay_tmp_{}.txt'.format(CMPD_LIMIT[i]))
		else:
			df_tmp = df_info.loc[(df_info[cmpd_num] >= CMPD_LIMIT[i]) & (df_info[cmpd_num] < CMPD_LIMIT[i - 1])]
			assay_class.append('assay_tmp_{}.txt'.format(CMPD_LIMIT[i]))
		
		assay_num_in_class.append(len(df_tmp))
		df_tmp.to_csv('{}/{}'.format(out_dir, assay_class[i]), header=False, sep='\t')
	

	df_tmp = df_info.loc[df_info[cmpd_num] < CMPD_LIMIT[-1]]
	assay_class.append('assay_tmp_st{}.txt'.format(CMPD_LIMIT[-1]))
	assay_num_in_class.append(len(df_tmp))
	df_tmp.to_csv('{}/{}'.format(out_dir, assay_class[-1]), header=False, sep='\t')
	
	#print('TEST: assay_num_in_class', assay_num_in_class)
	#print('TEST: assay_class', assay_class)
	
	assay_class = [a for i, a in enumerate(assay_class) if assay_num_in_class[i] > 0]

	mem = [m for i, m in enumerate(MDL_mem) if assay_num_in_class[i] > 0]
	
	if threads != 1:
		thd = [m for i, m in enumerate(threads) if assay_num_in_class[i] > 0]
	else:
		thd = [threads]
	assay_num_in_class = [a for a in assay_num_in_class if a > 0]

	
	return(assay_class, assay_num_in_class, mem, thd)

	

def index_master_table(input_file, out_file):
	# generate indices for each line of MasterTable or Z scaled files
	offset = 0
	#index starts from 0
	line_offset = [0]
	cmpds = []
	separator = ',' if input_file.endswith('.csv') else '\t'
	
	with open(input_file) as fid:	
		for line in fid:
			#print(i)
			offset += len(line)
			line24 = line[0 : 24]
			cmpds.append(line24.split(separator)[0])
			line_offset.append(offset)
	#remove the last index
	del line_offset[-1]
	
	df = pd.DataFrame([cmpds, line_offset]).T
	#df.to_csv(args.Output, index=False, header=False)
	df.set_index(df.columns[0], inplace=True)
	# save it to binary file
	joblib.dump(df, out_file)
	
def DataPreparation(pIC50file_ori, cmpdFile_ori, pIC50file, cmpdFile, PythonScriptFN):
	# validate the smiles structures.
	if os.path.isfile(pIC50file) and os.path.isfile(cmpdFile):
		print('\tCompounds already prepared')
	else:
		bash_line = get_bash(job_name='DataP', job_num=1, job_mem='32G', threads=1, py=True)
		bash_line += "python " + PythonScriptFN + " " + cmpdFile_ori + " " + pIC50file_ori + " " + cmpdFile + " " + pIC50file + "\n"

		bash_file = 'DataPrep.sh'
		with open(bash_file, 'w') as fid:
			fid.writelines(bash_line)
		
		process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free={}".format("32G"), bash_file, "/dev/null"], stdout=subprocess.PIPE)
		jobID = get_job_id(process)
		check_current_jobs([jobID])
		print('\tCompounds prepared')
		

def RF4pQSAR(resultsDir, pIC50file, PythonScriptFN, CMPD_LIMIT, RF_threads, RF_MDL_mem, PredvsExpRF, SummaryRF, RF_model_suffix, cmpd_num):
	
	# fraction of the training set compounds
	fraction2train = 0.75
	# the minimum number of compounds in one dataset
	MinMoleculesPerAssay = 50
	# the number of unique measurements in one dataset
	MOL_UNIQ_MIN = 5
	# TEST: 
	#MinMoleculesPerAssay = 45
	# 200 estimators in one random forest model
	trees = 200
	smiles = "smiles"
	CompoundIDcolumn = "compoundID"
	PIC50 = 'PIC50'
	AssayID = 'AssayID'
	info = "_info.txt"
	split = "_tmp_split.txt"
	predicted = "_pred.txt"
	SummaryFolder = '{}/{}'.format(resultsDir, 'Summary_pQSAR')
	
	if not os.path.isdir(SummaryFolder):
		os.mkdir(SummaryFolder)

	SummaryFile = '{}/{}'.format(SummaryFolder, 'SummaryRF.csv')
	AIDReference = '{}/{}'.format(SummaryFolder, 'AIDReference.csv')
	resultsDir = '{}/RF_Models'.format(resultsDir)

	if not os.path.isdir(resultsDir):
		os.mkdir(resultsDir)

	#If the summary file is there, all models were built
	if os.path.isfile(SummaryFile):
		print("\tNo models to do, summary file is there! Exiting...")
		sys.exit(1)

	if not os.path.isfile(AIDReference):			
		# AIDReference: a file contains assay ID and the number of compounds in associated assay ID
		df_all = pd.read_csv(pIC50file, header=0, index_col=None, sep=Tools.separator4fileName(pIC50file))
		colnames = list(df_all.columns)
		colnames[0 : 4] = [CompoundIDcolumn, AssayID, PIC50, smiles]
		df_all.columns = colnames
		df_info_all = pd.DataFrame(df_all.groupby(AssayID)[PIC50].count())
		df_info_all.columns = [cmpd_num]
		df_info_all['nunique'] = df_all.groupby(AssayID)[PIC50].nunique()
		df_info_all = df_info_all[(df_info_all[cmpd_num] >= MinMoleculesPerAssay) &  (df_info_all['nunique'] >= MOL_UNIQ_MIN)]
		df_info_all[cmpd_num].to_csv(AIDReference, index=True, header=False)
		del df_all
	
	else:
		df_info_all = pd.read_csv(AIDReference, index_col=0, header=None)
		df_info_all.columns = [cmpd_num]
	
	df_info_all.index = [str(i) for i in df_info_all.index]
	AssayIDs = list(df_info_all.index)
	completed = [f.split(info)[0] for f in os.listdir(resultsDir) if f.endswith(info)]
	AssayIDs = [item for item in AssayIDs if item not in completed]
	
	df_info = df_info_all.loc[AssayIDs]
	assay_class, assay_num_in_class, RF_MDL_mem, RF_threads = split_by_cmpd_number(df_info, CMPD_LIMIT, RF_MDL_mem, resultsDir, cmpd_num, RF_threads)
	
	print('\tStart building {} RF models'.format(sum(assay_num_in_class)))
	
	jobIDs = []
	for i, ac in enumerate(assay_class):
		## Submit jobs separately according to the number of compounds and allocate different memory.
		bash_line = get_bash(job_name='RFM', job_num=assay_num_in_class[i], job_mem=RF_MDL_mem[i], threads=RF_threads[i], py=True)
		bash_line += "SAMPLES_LIST='{}/{}'\n".format(resultsDir, ac)
		bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F'\\t' '{print $1}'`\n"
		bash_line += "echo -e '{}\\t{}\\t{}\\t{}'".format(CompoundIDcolumn, AssayID, PIC50, smiles) + ' > ' + resultsDir + '/' + '${ARGS}' + split + '\n'
		## assume pIC50FILE is tab separated
		bash_line += "awk -F'\\t' '$2==\"'${ARGS}'\"{print $1, $2, $3, $4}'" + ' ' + "OFS='\\t'" + ' ' + pIC50file + " >> " + resultsDir + "/$" + "{ARGS}" + split +"\n"

		#bash_line += "awk -F',' '$2==\"'${ARGS}'\"{print $1, $2, $3, $4}'" + ' ' + "OFS='\\t'" + ' ' + pIC50file + " >> " + resultsDir + "/$" + "{ARGS}" + split +"\n"
		bash_line += "python " + PythonScriptFN + " " + resultsDir + '/' + '${ARGS[0]}' + split + " " + resultsDir + '/' + '${ARGS[0]}' + "." + RF_model_suffix + " " + resultsDir + '/' + '${ARGS[0]}' + info + " " + resultsDir + '/' + '${ARGS[0]}' + predicted + " " + str(fraction2train) + " " + str(MinMoleculesPerAssay) + " " + CompoundIDcolumn + " " + smiles + " " + "TAB" + " " + str(trees) + " " + str(RF_threads[i]) + "\n"
		##TEST:
		#bash_line += 'rm ' + resultsDir + '/' + '${ARGS}' + split + '\n'
		bash_file = '{}/RF_MDL_{}_tmp.sh'.format(resultsDir, ac)
		with open(bash_file, 'w') as fid:
			fid.writelines(bash_line)
		
		process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free={}".format(RF_MDL_mem[i]), bash_file, "/dev/null"], stdout=subprocess.PIPE)
		jobID = get_job_id(process)
		jobIDs.append(jobID)
	
	check_current_jobs(jobIDs)

	# check for all completed models in the same directory
	completed = [f.split(info)[0] for f in os.listdir(resultsDir) if f.endswith(info)]
	failed = set(AssayIDs) - set(completed)

	if len(failed) > 0:
		print('\t{} RF models failed, due to some reasons'.format(len(failed)))
		AssayIDs = list(failed)
	else:
		print('\tAll models were built')

	#merge info and predicted vs. experimental
	subprocess.run(["awk 'NR==1 {header=$_} FNR==1 && NR!=1 { $_ ~ $header getline; } {print}' " + resultsDir + "/*" + info + " | tr '\\t' ',' > " + SummaryFile], shell=True)
	subprocess.run(["awk 'NR==1 {header=$_} FNR==1 && NR!=1 { $_ ~ $header getline; } {print}' " + resultsDir + "/*" + predicted + " | tr '\\t' ',' > " + resultsDir + "/" + PredvsExpRF], shell=True)

	#delete intermidiate files
	subprocess.run(["rm " + resultsDir + "/*" + info], shell=True)
	subprocess.run(["rm " + resultsDir + "/*" + predicted], shell=True)
	#subprocess.run(["rm " + resultsDir + "/*_tmp_*"], shell=True)
	print("\tFinished building RF models")


def RFpredict(resultsDir, file2predict, PythonScriptFN, modelDir, SPLIT_NUM, RF_PRED_mem, RF_PRED_threads, MasterTable):
	MasterTable_indices = 'MasterTable_indices'
	ResultExtension = "_master_table_tmp.csv"
	Job_name = "RFP"
	QC_name = 'RFP_QC'
	idx = False
	split_file = 'split_assay_tmp.txt'
	split_bash = 'split_bash_tmp.sh'
	suffix_Cmpd = '.Cmpd_list.csv'
	suffix_fps = '.morgan_fps.joblib'
	# adjust this variable to a larger for huge dataset
	chunk_size = 5_000
	CIDReference = '{}/{}/CIDReference.csv'.format(resultsDir, 'Summary_pQSAR')
	RF_dir = '{}/RF_Models'.format(resultsDir)
	resultsDir = '{}/RF_Predictions'.format(resultsDir)

	if not os.path.isdir(resultsDir):
		os.mkdir(resultsDir)
	
	QC_file = '{}/QC_RF_prediction_tmp.txt'.format(resultsDir)
	RFP_completed = '{}/RFP_completed_tmp.txt'.format(resultsDir)
	AID_MDL_Reference = '{}/AID_RF_MDL_Reference_tmp.csv'.format(resultsDir)

	if os.path.isfile('{}/{}'.format(resultsDir, MasterTable)):
		print("\tMasterTable is there! Exiting...")
		sys.exit(1)
	
	if not os.path.isfile(file2predict):
		print("\tNothing to predict... {} is not there.".format(file2predict))
		sys.exit(1)
	
	print("\tStart calculating Morgan fingerprint")
	smiles2fps(file2predict, resultsDir, suffix_Cmpd, suffix_fps, chunk_size)
	
	fps_list = ','.join(['{}/{}'.format(resultsDir, f) for f in os.listdir(resultsDir) if f.endswith(suffix_fps)])
	
	#generate Cmpd reference for final merging. 
	if not os.path.isfile(CIDReference):
		Cmpd_list = ['{}/{}'.format(resultsDir, f) for f in os.listdir(resultsDir) if f.endswith(suffix_Cmpd)]
		
		for i, Cmpd in enumerate(Cmpd_list):
			if i == 0:
				df_Cmpd = pd.read_csv(Cmpd, header=0, index_col=None)
			else:
				df_Cmpd = pd.concat([df_Cmpd, pd.read_csv(Cmpd, header=0, index_col=None)], axis=0, ignore_index=True)
		
		col_name = list(df_Cmpd.columns)
		df_Cmpd.sort_values(by=col_name, axis=0, inplace=True)
		df_Cmpd.to_csv('{}'.format(CIDReference), index=False, header=True)
	

	if not os.path.isfile('{}/{}'.format(resultsDir, split_file)):
		AID_mdl_list = [f.split('.')[0] for f in os.listdir(modelDir) if f.endswith('.RFmodel')]
		
		with open(AID_MDL_Reference, 'w') as fid:
			fid.writelines('\n'.join(AID_mdl_list) + '\n')
		
		bash_line = 'cd {}\n'.format(resultsDir)
		bash_line += "ASSAY='{}'\n".format(AID_MDL_Reference)
		bash_line += 'split -a3 -l{} $ASSAY tmp_\n'.format(SPLIT_NUM)
		## Three question marks. Corresponding to split -a3
		bash_line += 'ls tmp_??? > {}\n'.format(split_file)
		
		split_bash = '{}/{}'.format(resultsDir, split_bash)
		with open(split_bash, 'w') as fid:
			fid.writelines(bash_line)
		
		subprocess.run(["sh", split_bash], stdout=subprocess.PIPE)
	#this loop may restart after QCing merged CSV files
	QC_dict = {}
	
	print('\tStart prediction...')

	while True:
		## make RF prediction
		if os.path.getsize('{}/{}'.format(resultsDir, split_file)) == 0:
			break
		splits = pd.read_csv('{}/{}'.format(resultsDir, split_file), header=None, index_col=0)

		QC_job_num = len(splits)
		## For each loop, QC all existing predictions
		if os.path.isfile(QC_file):
			os.remove(QC_file)
		splits = list(splits.index)
		completed = [f.split(ResultExtension)[0] for f in os.listdir(resultsDir) if f.endswith(ResultExtension)]
		splits = list(set(splits) - set(completed))
		#TEST: print:
		
		with open('{}/{}'.format(resultsDir, split_file), 'w') as fid:
			fid.writelines('\n'.join(splits))
		if len(splits) > 0:
			print('\tStart predicting {} subsets...'.format(len(splits)))
			#get all model names
			
			bash_line = get_bash(job_name=Job_name, job_num=len(splits), job_mem=RF_PRED_mem, job_time=340000, py=True)

			bash_line += "SAMPLES_LIST='{}/{}'\n".format(resultsDir, split_file)
			bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F',' '{print $1}'`\n"
			bash_line += ' '.join(['python', PythonScriptFN, resultsDir + '/' + '${ARGS}', modelDir, 'RF', fps_list, resultsDir, str(idx), CIDReference]) + '\n'

			bash_file = '{}/{}_tmp.sh'.format(resultsDir, Job_name)
			with open(bash_file, 'w') as fid:
				fid.writelines(bash_line)
			
			process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free={}".format(RF_PRED_mem), bash_file, "/dev/null"], stdout=subprocess.PIPE)
			jobID = get_job_id(process)
			
			check_current_jobs([jobID])
			# TEST:
			#check_current_jobs([jobID], sec=60)
			
			# real time completed jobs
			completed2 = [f.split(ResultExtension)[0] for f in os.listdir(resultsDir) if f.endswith(ResultExtension)]
			with open(RFP_completed, 'w') as fid:
				fid.writelines('\n'.join(completed2) + '\n')
			
			if len(completed2) > 0:
				## QC the number of line for each successfully predicted file.
				print('\t Start quality checking...')
				bash_line = get_bash(job_name=QC_name, job_num=QC_job_num, job_mem='4G', job_time=340000, threads=1, py=False)
				bash_line += 'cd {}\n'.format(resultsDir)
				bash_line += "SAMPLES_LIST='{}'\n".format(RFP_completed)
				bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F',' '{print $1}'`\n"
				bash_line += "NUM=`wc -l " + resultsDir + "/" + '${ARGS}' + ResultExtension + '`\n'
				bash_line += 'echo $NUM >> {}\n'.format(QC_file)
				bash_file = '{}/{}_tmp.sh'.format(resultsDir, QC_name)
				
				with open(bash_file, 'w') as fid:
					fid.writelines(bash_line)
				
				process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free=4G", bash_file, "/dev/null"], stdout=subprocess.PIPE)
				jobID = get_job_id(process)
				
				check_current_jobs([jobID])
				# TEST:
				#check_current_jobs([jobID], sec=60)
				
				
				# Avoid empty lines or error lines (not started with number)
				with open(QC_file) as fid:
					for line in fid:
						if len(line.strip()) > 0:
							lf = line.split()
							try:
								num = int(lf[0])
								QC_dict[lf[1].split('/')[-1].split('_master_table_tmp')[0]] = num
							except:
								continue
							
				
				QC_NUM = QC_dict[max(QC_dict, key=QC_dict.get)]
				
				subprocess.run(["rm " + QC_file], shell=True)
				## list to delete
				failed = [k for k, v in QC_dict.items() if v < QC_NUM]
				
				## delete failed prediction and make a new split_file with failed names
				if len(failed) > 0:
					for f in failed:
						os.remove('{}/{}{}'.format(resultsDir, f, ResultExtension))
				
				print('\t\tQuality check done...')
		else:
			print('\t\tAll subsets were predicted')
			break
	
	## Merge predicted file
	print('\tStart merging predictions...')
	
	if not os.path.isfile('{}/{}'.format(resultsDir, MasterTable)):
		paste_pred(resultsDir, ResultExtension, CIDReference, MasterTable, Job_name='RFPST')
	print('\t\tMasterTable.csv generated...')
	## index the MasterTable.csv
	print('\tStart indexing MasterTable.csv...')
	
	if not os.path.isfile('{}/{}'.format(resultsDir, MasterTable_indices)):
		index_master_table('{}/{}'.format(resultsDir, MasterTable), '{}/{}'.format(resultsDir, MasterTable_indices))
	print('\t\tIndexing done...')

	subprocess.run(["rm " + resultsDir + "/*TMP*"], shell=True)
	subprocess.run(["rm " + RF_dir + "/*_tmp_*"], shell=True)


def PLS4pQSAR(resultsDir, MasterTableDir, PythonScriptFN, RF_expr_pred, CMPD_LIMIT, PLS_MDL_threads, PLS_MDL_mem, SummaryPLS, PLS_model_suffix, cmpd_num):
	# query profile by Cmpd ID on the fly
	split = "_Msplit.txt"
	info = "_PLSinfo.txt"
	predictedVexp = "_predictedPLS.txt"
	AssayColumnName = "PIC50"
	#placeHolder = "${ARGS[0]}"
	fractio2test = 0.75
	#headerName = "header.txt"
	smiles = "smiles"
	CompoundIDcolumn = "CompoundID"
	PIC50 = 'PIC50'
	AssayID = 'AssayID'
	
	SummaryFolder = '{}/{}'.format(resultsDir, 'Summary_pQSAR')
	SummaryFile = '{}/{}'.format(SummaryFolder, 'SummaryPLS.csv')

	#If the summary file is there, all models were built
	if os.path.isfile(SummaryFile):
		print("\tNo models to do, summary file is there! Exiting...")
		sys.exit(1)
	#########%%%%%%%%%%%%%%%%%%

	RF_mdl_dir = '{}/RF_Models'.format(resultsDir)
	resultsDir = '{}/PLS_Models'.format(resultsDir)

	if not os.path.isdir(resultsDir):
		os.mkdir(resultsDir)	

	AssayID_list = [f.split('.')[0] for f in os.listdir(RF_mdl_dir) if f.endswith('.RFmodel')]
		
	df_all = pd.read_csv('{}/{}'.format(SummaryFolder, 'SummaryRF.csv'), header=0, index_col=0)
	df_all.index = [str(i) for i in df_all.index]
	df_all = df_all[['Count(molecules)']]
	df_all.columns = [cmpd_num]
		
	df_info_all = df_all.loc[set(df_all.index).intersection(set(AssayID_list))].copy()
	##################%%%%%%%%%%%%%%%%%%%%%

	AssayIDs = list(df_info_all.index)
	completed = [f.split(info)[0] for f in os.listdir(resultsDir) if f.endswith(info)]
	AssayIDs = [item for item in AssayIDs if item not in completed]
	

	while True:
		df_info = df_info_all.loc[AssayIDs]
		assay_class, assay_num_in_class, PLS_MDL_mem, PLS_threads = split_by_cmpd_number(df_info, CMPD_LIMIT, PLS_MDL_mem, resultsDir, cmpd_num, PLS_MDL_threads)
		len_list = [len(assay_class), len(assay_num_in_class), len(PLS_MDL_mem), len(PLS_threads)]

		print('\tStart building {} PLS models'.format(sum(assay_num_in_class)))
		jobIDs = []
		for i, ac in enumerate(assay_class):
			## Submit jobs separately according to the number of compounds. Allocate different memory.
			bash_line = get_bash(job_name='PLSM', job_num=assay_num_in_class[i], job_mem=PLS_MDL_mem[i], threads=PLS_threads[i], py=True)
			bash_line += "SAMPLES_LIST='{}/{}'\n".format(resultsDir, ac)
			bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F'\\t' '{print $1}'`\n"
			bash_line += "echo -e '{}\\t{}\\t{}\\t{}\\t{}'".format(CompoundIDcolumn, AssayID, PIC50, smiles, 'Set') + ' > ' + resultsDir + '/' + '${ARGS}' + split + '\n'
			bash_line += "awk -F',' '$2==\"'${ARGS}'\"{print $1, $2, $3, $4, $5}'" + ' ' + "OFS='\\t'" + ' ' + RF_expr_pred + " >> " + resultsDir + "/$" + "{ARGS}" + split +"\n"
			#bash_line += "python {} {} ${ARGS[0]}{} {} ${ARGS[0]} ${ARGS[0]}.{} ${ARGS[0]}{} ${ARGS[0]}{} {}".format(PythonScriptFN, MasterTableDir, split, AssayColumnName, PLS_model_suffix, info, predictedVexp, str(fractio2test))
			bash_line += 'python ' + PythonScriptFN + ' ' +  MasterTableDir + ' ' + resultsDir + '/' + '${ARGS}' + split + ' ' + AssayColumnName + ' ${ARGS} ' + resultsDir + '/' + '${ARGS}' + '.' + PLS_model_suffix + ' ' + resultsDir + '/' + '${ARGS}' + info + ' ' + resultsDir + '/' + '${ARGS}' + predictedVexp + ' ' + str(fractio2test) + '\n'
			##TEST:
			bash_line += 'rm ' + resultsDir + '/' + '${ARGS}' + split + '\n'	
			bash_file = '{}/PLS_MDL_{}_tmp.sh'.format(resultsDir, ac)
			with open(bash_file, 'w') as fid:
				fid.writelines(bash_line)
			
			process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free={}".format(PLS_MDL_mem[i]), bash_file, "/dev/null"], stdout=subprocess.PIPE)
			jobID = get_job_id(process)
			jobIDs.append(jobID)
		
		check_current_jobs(jobIDs)

		#check for all completed models in the same directory
		completed = [f.split(info)[0] for f in os.listdir(resultsDir) if f.endswith(info)]
		failed = set(AssayIDs) - set(completed)

		if len(failed) > 0:
			print('\t{} PLS models failed, will try again'.format(len(failed)))
			AssayIDs = list(failed)
		else:
			print('\t\tAll models were built')
			break
	
	#2. Merge all files and delete temporary files
	print('\tMerging summary files')
	MergedPredName = "PredvsExpPLS.csv"
	command = "awk 'NR==1 {header=$_} FNR==1 && NR!=1 { $_ ~ $header getline; } {print}' " + resultsDir + '/' + "*" + info + " | tr '\\t' ','" + " > " + SummaryFile
	subprocess.run([command], shell=True)
	command = "awk 'NR==1 {header=$_} FNR==1 && NR!=1 { $_ ~ $header getline; } {print}' " + resultsDir + '/' + "*" + predictedVexp +  " | tr '\t' ','" + " > " + resultsDir + '/' + MergedPredName
	subprocess.run([command], shell=True)

	subprocess.run(["rm " + resultsDir + '/' + "*" + info], shell=True)
	subprocess.run(["rm " + resultsDir + '/' + "*" + predictedVexp], shell=True)
	subprocess.run(["rm " + resultsDir + '/' + "*tmp*"], shell=True)
	print('\t\tDone...')


def PLSpredict(resultsDir, MasterTableFN, PythonScriptFN, modelDir, CIDReference, SPLIT_NUM, PLS_PRED_mem, ZscaledAllAssaysCID, ZscaledAllAssaysAID, ZscalingStats):
	
	ZscaledAllAssaysCID = 'ZscaledAllAssaysCID.csv'
	ZscaledAllAssaysAID = 'ZscaledAllAssaysAID.csv'
	ZscaledAllAssaysCID_indices = 'ZscaledAllAssaysCID_indices'
	ZscaledAllAssaysAID_indices = 'ZscaledAllAssaysAID_indices'	
	ZscalingStats = 'ZscalingStats.csv'
	ZscalingExtension = '_stat_table_tmp.csv'
	CmpdExtension = "_master_table_tmp_CID.csv"
	AIDExtension = "_master_table_tmp_AID.csv"
	
	Job_name = "PLSP"
	Cmpd_QC = 'Cmpd_QC'

	idx = False
	split_file = 'split_assay_tmp.txt'
	split_bash = 'split_bash_tmp.sh'
	PLS_mdl_dir = '{}/PLS_Models'.format(resultsDir)
	resultsDir = '{}/PLS_Predictions'.format(resultsDir)

	if not os.path.isdir(resultsDir):
		os.mkdir(resultsDir)
	
	AID_PLS_MDL_ref = [AID.split('.')[0] for AID in os.listdir(PLS_mdl_dir) if AID.endswith('.PLSmodel')]
	
	AID_MDL_Reference = '{}/AID_PLS_MDL_Reference_tmp.csv'.format(resultsDir)
	with open (AID_MDL_Reference, 'w') as fid:
		fid.writelines(['\n'.join(AID_PLS_MDL_ref)])

	if os.path.isfile('{}/{}'.format(resultsDir, ZscaledAllAssaysCID)) and os.path.isfile('{}/{}'.format(resultsDir, ZscaledAllAssaysAID)):
		print('\tAll subsets were predicted')
	else:
		#Get available PLS models
		# resume from last abrupt point  
		if not os.path.isfile('{}/{}'.format(resultsDir, split_file)):
			AID_mdl_list = [f.split('.')[0] for f in os.listdir(modelDir) if f.endswith('.PLSmodel')]
			
			with open(AID_MDL_Reference, 'w') as fid:
				fid.writelines('\n'.join(AID_mdl_list) + '\n')
			
			
			bash_line = 'cd {}\n'.format(resultsDir)
			bash_line += "ASSAY='{}'\n".format(AID_MDL_Reference)
			bash_line += 'split -a3 -l{} $ASSAY tmp_\n'.format(SPLIT_NUM)
			## Three question marks. Corresponding to split -a3
			bash_line += 'ls tmp_??? > {}\n'.format(split_file)
			
			split_bash = '{}/{}'.format(resultsDir, split_bash)
			with open(split_bash, 'w') as fid:
				fid.writelines(bash_line)
			
			subprocess.run(["sh", split_bash], stdout=subprocess.PIPE)
		
		QC_dict = {}
		
		while True:
			## make PLS prediction
			try:
				splits = pd.read_csv('{}/{}'.format(resultsDir, split_file), header=None, index_col=0)
				splits = list(splits.index)
				#print('splits', splits)
				completed = [f.split(CmpdExtension)[0] for f in os.listdir(resultsDir) if f.endswith(CmpdExtension)]
				splits = list(set(splits) - set(completed))
				# save the splits that have not been calculated.
				with open('{}/{}'.format(resultsDir, split_file), 'w') as fid:
					fid.writelines('\n'.join(splits))
			except:
				splits = []
			
			if len(splits) > 0:
				print('\tStart predicting {} subsets...'.format(len(splits)))
				#get all model names
				
				bash_line = get_bash(job_name=Job_name, job_num=len(splits), job_mem=PLS_PRED_mem, job_time=340000, py=True)

				bash_line += "SAMPLES_LIST='{}/{}'\n".format(resultsDir, split_file)
				bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F',' '{print $1}'`\n"
				bash_line += ' '.join(['python', PythonScriptFN, resultsDir + '/' + '${ARGS}', modelDir, 'PLS', MasterTableFN, resultsDir, str(idx), CIDReference]) + '\n'
				bash_file = '{}/{}_tmp.sh'.format(resultsDir, Job_name)
				with open(bash_file, 'w') as fid:
					fid.writelines(bash_line)
				
				process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free={}".format(PLS_PRED_mem), bash_file, "/dev/null"], stdout=subprocess.PIPE)
				jobID = get_job_id(process)
				
				check_current_jobs([jobID])
				
				## QC the number of line for each predicted file.
				## QC the long thin (Cmpd_QC) and short wide (AID_QC) Z scaled files 
				
				print('\tStart Quality checking ...')
				
				#for QC in [Cmpd_QC, AID_QC]:
				for QC in [Cmpd_QC]:
					bash_line = get_bash(job_name=QC, job_num=len(splits), job_mem='4G', job_time=340000, threads=1, py=False)
					bash_line += 'cd {}\n'.format(resultsDir)
					bash_line += "SAMPLES_LIST='{}/{}'\n".format(resultsDir, split_file)
					bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F',' '{print $1}'`\n"
					bash_line += "NUM=`wc -l " + resultsDir + "/" + '${ARGS}' + CmpdExtension + '`\n'
					bash_line += 'echo $NUM >> {}/{}_tmp.txt'.format(resultsDir, QC)
					
					bash_file = '{}/{}_tmp.sh'.format(resultsDir, QC)
					with open(bash_file, 'w') as fid:
						fid.writelines(bash_line)
					
					process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free=4G", bash_file, "/dev/null"], stdout=subprocess.PIPE)
					jobID = get_job_id(process)
					
					check_current_jobs([jobID])
			
					with open('{}/{}_tmp.txt'.format(resultsDir, QC)) as fid:
						for line in fid:
							if len(line.strip()) > 0:
								lf = line.split()
								try:
									num = int(lf[0])
									QC_dict[lf[1].split('/')[-1].split('_master_table_tmp')[0]] = num
								except:
									continue
					
						if QC == Cmpd_QC:
							QC_NUM = QC_dict[max(QC_dict, key=QC_dict.get)]
							## list to delete
							Cmpd_failed = [k for k, v in QC_dict.items() if v < QC_NUM]
						
				subprocess.run(["rm " + resultsDir + '/' + "{}_tmp.txt".format(QC)], shell=True)
				failed = list(set(Cmpd_failed))
				## delete failed prediction and make a new split_file with failed names
				if len(failed) > 0:
					for f in failed:
						os.remove('{}/{}{}'.format(resultsDir, f, CmpdExtension))
						os.remove('{}/{}{}'.format(resultsDir, f, AIDExtension))
				print('\t\tQuality check done...')
			else:
				break
		print('\tStart merging Cmpd predictions...')
		
		if not os.path.isfile('{}/{}'.format(resultsDir, ZscaledAllAssaysCID)):
			paste_pred(resultsDir, CmpdExtension, CIDReference, ZscaledAllAssaysCID, Job_name='NVPPST')
		
		print('\tStart merging AID predictions...')
		if not os.path.isfile('{}/{}'.format(resultsDir, ZscaledAllAssaysAID)):
			paste_pred2(resultsDir, AIDExtension, ZscaledAllAssaysAID, Job_name='AIDPST')
			
		print('\tStart merging Z scaling stats...')
		if not os.path.isfile('{}/{}'.format(resultsDir, ZscalingStats)):
			paste_pred2(resultsDir, ZscalingExtension, ZscalingStats, Job_name='StatPST')
		
		print('\t\tMasterTable generated...')
		
	print('\tStart indexing MasterTable...')
	
	if not os.path.isfile('{}/{}'.format(resultsDir, ZscaledAllAssaysCID_indices)):
		index_master_table('{}/{}'.format(resultsDir, ZscaledAllAssaysCID), '{}/{}'.format(resultsDir, ZscaledAllAssaysCID_indices))
	
	if not os.path.isfile('{}/{}'.format(resultsDir, ZscaledAllAssaysAID_indices)):
		index_master_table('{}/{}'.format(resultsDir, ZscaledAllAssaysAID), '{}/{}'.format(resultsDir, ZscaledAllAssaysAID_indices))
	print('\t\tIndexing done...')
	subprocess.run(["rm " + resultsDir + "/*tmp*"], shell=True)
	
def IndexTable(MasterTableFN, IndexFN, PythonScriptFN):
	# Generate binary index (joblib)
	Directory, _ = os.path.split(MasterTableFN)
	bash_line = get_bash(job_name='Index', job_mem='16G', job_time=340000, threads=1, py=True)
	bash_line += 'cd {}\n'.format(Directory)
	bash_line += 'python {} -i {} -o {}\n'.format(PythonScriptFN, MasterTableFN, IndexFN)

	bash_file = '{}/index_tmp.sh'.format(Directory)
	
	process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free=16G", bash_file, "/dev/null"], stdout=subprocess.PIPE)
	
	jobID = get_job_id(process)
	check_current_jobs([jobID])

def index_marks(nrows, chunk_size):
    return range(1 * chunk_size, (nrows // chunk_size + 1) * chunk_size, chunk_size)

def split_Cmpd(dfm, chunk_size):
    indices = index_marks(dfm.shape[0], chunk_size)
    return np.split(dfm, indices)


def py4fps(py4fps_tmp):
	py = 'import os\n'
	py += 'import sys\n'
	py += 'import pandas as pd\n'
	py += 'import joblib\n'
	py += 'from rdkit.Chem import AllChem, PandasTools\n'	

	py += "smiles = 'smiles'\n"
	py += "MOLECULE = 'MOLECULE'\n"
	py += "FP = 'FP'\n"	
	py += "separator = ',' if sys.argv[1].endswith('.csv') else '\t'\n"
	py += 'df_cmpd = pd.read_csv(sys.argv[1], index_col=0, sep=separator)\n'
	py += "df_cmpd.columns = ['smiles']\n"
	py += 'PandasTools.AddMoleculeColumnToFrame(df_cmpd, smiles, MOLECULE)\n'
	py += 'df_cmpd = df_cmpd.loc[df_cmpd[MOLECULE].notnull()]\n'
	# example for morgan radius 2, 1024-bit fingerprint
	py += 'df_cmpd[FP] = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in df_cmpd[MOLECULE]]\n'	
	py += 'df_cmpd.drop([MOLECULE, smiles], axis=1, inplace=True)\n'
	py += 'df_cmpd.dropna(axis=0, subset=[FP], inplace=True)\n'	
	py += "df_cmpd = pd.DataFrame([list(x) for x in df_cmpd[FP]], columns=list(range(1024)), index=df_cmpd.index, dtype='int8')\n"
	py += 'joblib.dump(df_cmpd, sys.argv[1] + sys.argv[3])\n'
	py += "df_cmpd.index.name = 'compoundID'\n"
	py += 'df_cmpd = pd.DataFrame(df_cmpd.index)\n'
	py += "df_cmpd.to_csv(sys.argv[1] + sys.argv[2], index=False)\n"
	
	with open(py4fps_tmp, 'w') as fid:
		fid.writelines(py)


def smiles2fps(input_file, out_dir, suffix_Cmpd, suffix_fps, chunk_size):
	# Generate fps
	py4fps_tmp = '{}/py4fps_tmp.py'.format(out_dir)
	chunks_list_tmp = '{}/chunks_list_TMP.txt'.format(out_dir)
	separator = Tools.separator4fileName(input_file)
	df_cmpd = pd.read_csv(input_file, header=0, index_col=0, sep=separator)
	df_cmpd.dropna(axis=0, inplace=True)
	chunks = split_Cmpd(df_cmpd, chunk_size)
	chunks_list = []
	
	for i, chk in enumerate(chunks):
		tmp_name = '{}/chunk_{}_Cmpd_TMP.csv'.format(out_dir, i)
		chk.to_csv(tmp_name , header=True, index=True)
		chunks_list.append(tmp_name)
			
	completed = ['{}/{}'.format(out_dir, f.split(suffix_fps)[0]) for f in os.listdir(out_dir) if f.endswith(suffix_fps)]
	#TEST:
	#print('CHUNK completed', completed)
	chunks_list = list(set(chunks_list) - set(completed))
	#TEST:
	#print('CHUNK tmp_name', chunks_list)
	if len(chunks_list) > 0:
	
		with open(chunks_list_tmp, 'w') as fid:
			fid.writelines('\n'.join(chunks_list) + '\n')

		py4fps(py4fps_tmp)
		
		bash_line = get_bash(job_name='FPS', job_num=len(chunks_list), job_mem='16G', job_time=340000, threads=1, py=True)
		bash_line += 'cd {}\n'.format(out_dir)
		bash_line += "SAMPLES_LIST='{}'\n".format(chunks_list_tmp)
		bash_line += "ARGS=`sed -n \"${SGE_TASK_ID}p\" < $SAMPLES_LIST | awk -F'\\t' '{print $1}'`\n"
		bash_line += "python " + py4fps_tmp + ' ' + '${ARGS[0]}' + ' ' + suffix_Cmpd + ' ' + suffix_fps + ' ' + out_dir + '\n'

		bash_file = '{}/bash4fps_tmp.sh'.format(out_dir)
		with open(bash_file, 'w') as fid:
			fid.writelines(bash_line)
		process = subprocess.run(["qsub", "-l", "h_rt=340000" + ",m_mem_free=16G", bash_file, "/dev/null"], stdout=subprocess.PIPE)
		jobID = get_job_id(process)
		
		check_current_jobs([jobID])


def pQSARsummary(RFsummary, PLSsummary, ZscaleStats, AssaySummary, PIC50FILE, resultFile):
	#
	general_info = '{}/SummaryPQSAR_all.csv'.format(resultFile)
	simplified_info = '{}/SummaryPQSAR.csv'.format(resultFile)

	if os.path.isfile(simplified_info):
		print('\tSummary file is there! Exiting...')
		sys.exit(1)


	df_rf = pd.read_csv(RFsummary, index_col=0, header=0, sep=',')
	df_rf.dropna(how='all', axis=1, inplace=True)
	
	df_pls = pd.read_csv(PLSsummary, index_col=0, header=0, sep=',')
	df_stats = pd.read_csv(ZscaleStats, index_col=0, header=0, sep=',')
	df_rf = df_rf.merge(df_pls, right_index=True, left_index=True, how='right')
	df_rf = df_rf.merge(df_stats, right_index=True, left_index=True, how='left')

	# first column is CDS_ID and second ASSAY_VER_ID
	df_assay = pd.read_csv(AssaySummary, index_col=None, header=0, sep=',')
	df_assay = df_assay.groupby('assay_chembl_id').nth(0)
	
	df_data = pd.read_csv(PIC50FILE, index_col=None, header=0, sep='\t' if PIC50FILE.endswith('.txt') else ',')
	
	df_group = df_data[['AID', 'standard_units', 'standard_type', 'transform']].groupby('AID').nth(0)
	#df_group = df_data[['CDS_ID', 'standard_units', 'standard_type', 'transform']].groupby('CDS_ID').nth(0)
	
	df_assay = df_assay.merge(df_group, how='left', right_index=True, left_index=True)

	df_rf = df_rf.merge(df_assay, right_index=True, left_index=True, how='left')
	df_rf.index.name = 'AssayID'
	df_rf = df_rf[~df_rf.index.duplicated(keep='first')]
	df_rf.fillna('--', inplace=True)
	df_rf.to_csv(general_info, index = True, header=True,  float_format = "%.2f", sep=Tools.separator4fileName(general_info))

	cols = ["R^2_RF(ext)", "Count(molecules)", "R^2_RF(random)", "Actives(>=6)", "pIC50_min", "pIC50_max", "STD(pIC50)", "R2ext", "NV", "Q2orig", "R^2fit", "Compounds#", "Threshold", "Columns#", "stdev_pred", "mean_pred", "max_pred", "min_pred", "standard_units", "standard_type", "transform", "target_chembl_id", "description", "assay_type", "site_name", "target_type", "pref_name", "organism", "assay_organism", "assay_strain", "assay_tissue", "assay_cell_type", "confidence_score"]

	df_rf[cols].to_csv(simplified_info, index = True, header=True,  float_format = "%.2f", sep=Tools.separator4fileName(general_info))		


def main():	
	
	parser = argparse.ArgumentParser(description='Program for parsing a file wiht Z-scaled pQSAR predictions:', epilog="")

	# add long and short arguments
	parser.add_argument("--input", "-i", help="File with compound IDs. Only the first column will be used and interpreted as IDs")
	parser.add_argument("--output", "-o", help="File name for outputing results. The Z-scaled matrix will be trasposed")
	parser.add_argument("--file1", "-f", help="First additional file")
	parser.add_argument("--file2", "-f2", help="additional file")
	parser.add_argument("--file3", "-f3", help="additional file")
	parser.add_argument("--file4", "-f4", help="additional file")
	parser.add_argument("--job", "-j", help="Execute an SQL statement in the input file and save results to the output file")
	parser.add_argument("--CIDReference", "-cpd", help="compound ID list")
	parser.add_argument("--dir", "-d", help="Directory (normally working dir)")
	parser.add_argument("--pythonScript", "-p", help="Path and file name of a (python) script to be called from this application(if no path, the script will be in the same folder as this app)")
	
	########################################################
	if len(sys.argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	########################################################
	CIDReference = 'CIDReference.csv'
	ZscalingStats = 'ZscalingStats.csv'
	PredvsExpRF = 'PredvsExp.csv'
	SummaryRF = 'SummaryRF.csv'

	SummaryPLS = 'SummaryPLS.csv'
	SummaryRFR = 'SummaryRFR.csv'

	MasterTable = 'MasterTable.csv'
	#MasterTable_indices = 'MasterTable_indices'
	ZscaledAllAssaysCID = 'ZscaledAllAssaysCID.csv'
	ZscaledAllAssaysAID = 'ZscaledAllAssaysAID.csv'

	## Suffix
	RF_model_suffix = 'RFmodel'
	PLS_model_suffix = 'PLSmodel'
	cmpd_num = '#cmpd'
	
	## Parameters
	CMPD_LIMIT = [60_000, 20_000, 10_000, 5000, 2000]
	RF_MDL_threads = [4, 4, 4, 2, 1, 1]
	RF_MDL_mem = ['12G', '6G', '4G', '3G', '4G', '2G']
	# 8G is absolutely enough for RF prediction of Cmpd (actual observation is less than 4G). 
	RF_PRED_mem = '12G'
	RF_PRED_threads = 1
	PLS_MDL_threads = [4, 4, 2, 1, 1, 1]
	PLS_MDL_mem = ['16G', '8G', '8G', '8G', '4G', '3G']
	PLS_PRED_mem = '36G'
	SPLIT_NUM = 25
	# read arguments from the command line
	args = parser.parse_args()
	# check for --input
	if not args.input:
		print("Error: no input file provided")
		sys.exit(1)

	# check for --output
	if not args.output:
		print("Error: no output file to save results")
		sys.exit(1)

	# check for --job
	if not args.job:
		print("Error: no job is specified")
		sys.exit(1)
	if args.job == 'DataPrep':
		print('STEP 0: Data preparation...')
		pIC50file_ori = args.input
		cmpdFile_ori = args.file1
		pIC50file = args.file2
		cmpdFile = args.output
		PythonScriptFN = args.pythonScript
		DataPreparation(pIC50file_ori, cmpdFile_ori, pIC50file, cmpdFile, PythonScriptFN)


	if args.job == 'RF4pQSAR':
		print('STEP 1: Building Random Forest Regression Models...')
		resultsDir = args.output
		pIC50file = args.input
		PythonScriptFN = args.pythonScript
		
		RF4pQSAR(resultsDir, pIC50file, PythonScriptFN, CMPD_LIMIT, RF_MDL_threads, RF_MDL_mem, PredvsExpRF, SummaryRF, RF_model_suffix, cmpd_num)

	if args.job == 'RFpredict':
		print('STEP 2: RFR Model Prediction and Indexing...')
		# check for --input
		if not args.dir:
			print("\tError: no model directory path is specified!")
			sys.exit(1)

		modelDir = args.dir
		resultsDir = args.output
		file2predict = args.input
		PythonScriptFN = args.pythonScript

		RFpredict(resultsDir, file2predict, PythonScriptFN, modelDir, SPLIT_NUM, RF_PRED_mem, RF_PRED_threads, MasterTable)

	if args.job == 'PLS4pQSAR':
		print('STEP 3: Building PLS Regression Models...')
		resultsDir = args.output
		RF_expr_pred = args.input
		MasterTableDir = args.dir
		PythonScriptFN = args.pythonScript
		PLS4pQSAR(resultsDir, MasterTableDir, PythonScriptFN, RF_expr_pred, CMPD_LIMIT, PLS_MDL_threads, PLS_MDL_mem, SummaryPLS, PLS_model_suffix, cmpd_num)
		
	if args.job == 'PLSpredict':
		print('STEP 4: PLS Model Prediction and Indexing...')
		# check for --input
		modelDir = args.dir
		resultsDir = args.output
		MasterTableFN = args.input
		CIDReference = args.CIDReference
		PythonScriptFN = args.pythonScript

		PLSpredict(resultsDir, MasterTableFN, PythonScriptFN, modelDir, CIDReference, SPLIT_NUM, PLS_PRED_mem, ZscaledAllAssaysCID, ZscaledAllAssaysAID, ZscalingStats)
	
	#read file with compound IDs to search
	#generates CSV summary file
	if args.job == 'pQSARsummary':
		print('STEP 5: Summarizing...')
		# 
		# -i ${MDLADDR}/Summary_pQSAR/SummaryRF.csv -f ${MDLADDR}/Summary_pQSAR/SummaryPLS.csv -f2 ${MDLADDR}/PLS_Predictions/ZscalingStats.csv
		RFsummary = args.input
		PLSsummary = args.file1
		ZscaleStats = args.file2
		AssaySummary = args.file3
		PIC50FILE = args.file4
		resultFile = args.output
		pQSARsummary(RFsummary, PLSsummary, ZscaleStats, AssaySummary, PIC50FILE, resultFile)
	
if __name__ == "__main__":
	import os
	import sys
	import time
	import re
	import argparse
	import subprocess
	import numpy as np
	import pandas as pd
	import joblib
	import CommonTools as Tools
	main()
	
