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

def partition (list_in, n):
    #random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

def query_sql(activities, assay_chembl_id_new, preference_first, preference_second, act_sfile, output_lock):
    filelock = output_lock

    for i, assay in enumerate(assay_chembl_id_new):
        df_assay = activities[activities['assay_chembl_id'] == assay].copy()

        """
        First to capture all the assays with pchembl values: len(df_assay.pchembl.unique()) > 1
        After capturing assays with pchembl values. Then apply standard procedure for other assays with endpoint like
        EC50, CC50: len(set(df_assay['pchembl'])) <= 1.
        """

        if len(df_assay.pchembl.unique()) > 1:
            # df_assay.pchembl.unique() will be the unqiue pchembl values. len(df_assay.pchembl.unique())
            df_assay['PIC50'] = df_assay['pchembl']
            df_assay['transform'] = 'pchembl'
        else:
            # df_assay.pchembl will be a series of Nan if one assay has no pchembl values and df_assay.pchembl.unique() will be array([nan]) with length of 1. 
            assay_type = list(df_assay['standard_type'].unique())
            FLAG = True

            for tp in preference_first:
                if tp in assay_type:
                    break
            else:
                FLAG = False
            
            if FLAG == True:
                df_assay = df_assay[df_assay['standard_type'] == tp]
                df_assay = df_assay[df_assay['standard_relation'] == '=']
                if len(df_assay) >= 50:
                    # reverse the sign
                    if sum(df_assay['PIC50'] > 0) < 10:
                        df_assay['PIC50'] = - df_assay['PIC50']
                        # save df_assay
                else:
                    continue

            if FLAG == False:
                for tp in preference_second:
                    if tp in assay_type:
                        break
                else:
                    continue
                df_assay = df_assay[df_assay['standard_type'] == tp]
                
                if tp in ['IC50 ratio', 'ID50 ratio', 'Ratio EC50', 'Ratio_Papp']:
                    df_assay['PIC50'] = np.log10(df_assay['PIC50'])
                    df_assay.loc[df_assay['standard_relation'].isin(['>', '>=']), 'PIC50'] = df_assay.loc[df_assay['standard_relation'].isin(['>', '>=']), 'PIC50'].values + 1
                    df_assay.loc[df_assay['standard_relation'].isin(['<', '<=']), 'PIC50'] = df_assay.loc[df_assay['standard_relation'].isin(['<', '<=']), 'PIC50'].values - 1
                else:
                    units = list(df_assay['standard_units'].unique())
                    
                    if len(units) > 0:

                        if len(units) > 1:
                            units_dict = {u: len(df_assay[df_assay['standard_units'] == u]) for u in units}
                            unit = max(units_dict, key=units_dict.get)
                        else:
                            unit = units[0]

                        if unit in ['mg.kg-1', 'mg l-1', 'ug ml-1', 'ug.mL-1']:
                            mw = df_assay['mw']
                            transform = 1/(mw * 1000)
                            df_assay['PIC50'] = -np.log10(df_assay['PIC50'] * transform)
                            df_assay['transform'] = '`-log10<-{}'.format(unit)
                        elif unit == 'mM':
                            transform = 1/1000
                            df_assay['PIC50'] = -np.log10(df_assay['PIC50'] * transform)
                            df_assay['transform'] = '`-log10<-{}'.format(unit)
                        elif unit == 'nM':
                            transform = 1e-9
                            df_assay = df_assay[df_assay['PIC50'] > 0]
                            df_assay['PIC50'] = -np.log10(df_assay['PIC50'] * transform)
                            df_assay['transform'] = '`-log10<-{}'.format(unit)
                        elif unit in ['um', 'uM', 'uM.hr/L', 'umol.kg-1']:
                            transform = 1e-6
                            df_assay['PIC50'] = -np.log10(df_assay['PIC50'] * transform)
                            df_assay['transform'] = '`-log10<-{}'.format(unit)
                        
                        else:
                            continue
            
                    df_assay.loc[df_assay['standard_relation'].isin(['>', '>=']), 'PIC50'] = df_assay.loc[df_assay['standard_relation'].isin(['>', '>=']), 'PIC50'].values + 1
                    df_assay.loc[df_assay['standard_relation'].isin(['<', '<=']), 'PIC50'] = df_assay.loc[df_assay['standard_relation'].isin(['<', '<=']), 'PIC50'].values - 1

            
        df_assay.dropna(axis=0, subset=['PIC50', 'canonical_smiles'], inplace=True)
        df_assay.drop('pchembl', axis=1, inplace=True)
        stdev = np.std(df_assay['PIC50'])
        if len(df_assay) >= 50 and stdev >= 0.50 and stdev <= 10:
            if os.path.isfile(act_sfile):
                headers = False
            else:
                headers = True
            with filelock:
                df_assay.round(3).to_csv(act_sfile, mode='a', header=headers, index=False, sep='\t')
def main():

    dbase = sys.argv[1]
    addr = sys.argv[2]
    version = sys.argv[3]
    assay_file = sys.argv[4]
    threads = int(sys.argv[5])

    conn = sqlite3.connect(dbase)
    

    assay_sfile = '{}/{}_AllAssayInfo.csv'.format(addr, version)
    cmpd_sfile = '{}/{}_AllAssayCmpd.csv'.format(addr, version)
    act_sfile = '{}/{}_AllAssayData.txt'.format(addr, version)

    preference_first = [ "pEC50", "pI50", "pIC50", "pKb", "pKD", "pKi", "pKm", "pA2(app)", "pRA1","pRA2","log(10^6/IC50)", "log(1/C)", "Log 1/C", "log(1/I50)", "Log 1/Ki app", \
        "Log 1/Km", "Log 1/ratio", "log(activity)", "logEC50", "Log EC50", "Log I50", "log(IA)", "logIC50", "Log IC50", "-Log KB", "logKi", "Log Ki", "Log koff", "Log kon", "logPF", \
        "logRA", "log(RBA)", "Log 1/Km", "Km", "Log kon", "Log koff", "log(1/I50)", "log(1/C)", "Log 1/C", "Log 1/Ki app", "log(10^6/IC50)", "Relative binding affinity", "Thermal melting change",\
        "Delta G", "-Delta G", "Change in Cl- current", "Cl - current change", "Activity_index"]
    preference_second = ["k_off", "Ki", "Kd", "KA", "Ke", "Ke(app)", "AC50", "IC50", "EC50", "ED50", "CC50", "MIC", "Potency", "Potency_Clk2_uM", "Potency_Clk4_uM", \
        "Potency_Dyrk1a_uM", "Potency_Dyrk1b_uM","IC20", "IC90", "IC50 ratio", "ID50 ratio", \
        "Residual Activity", "Imax", "Beta2 duration", "Emax", "Inhibition", "Kd apparent", "KI_MICROM", "%max", \
        "max activation", "Papp", "Ratio EC50", "Ratio_Papp", "Response", "SS"]


    if 'v24' in version:
        assay_query = "select A.chembl_id as assay_chembl_id, A.assay_id, T.chembl_id as target_chembl_id, A.description, A.assay_type, B.site_name, T.target_type, T.pref_name, T.organism,  A.assay_organism, A.assay_strain, A.assay_tissue, A.assay_cell_type, A.confidence_score from assays A left join binding_sites B on A.tid == B.tid left join target_dictionary T on A.tid = T.tid"
        #cmpd_query = "select md.chembl_id, cs.canonical_smiles, cs.standard_inchi_key from compound_structures cs, molecule_dictionary md where cs.molregno == md.molregno;"
        cmpd_query = "select md.chembl_id as cmpd_chembl_id, cs.canonical_smiles, cp.mw_freebase as mw from compound_structures cs, molecule_dictionary md, compound_properties cp where cs.molregno == md.molregno and cs.molregno == cp.molregno;"
        #activity_query = "select md.chembl_id as cmpd_chembl_id, AY.chembl_id as assay_chembl_id,  AC.standard_value, AC.standard_relation, AC.standard_units, AC.standard_type, AC.data_validity_comment from activities AC, molecule_dictionary md, assays AY where AC.assay_id == AY.assay_id and AC.molregno == md.molregno and AC.data_validity_comment IS NOT 'Outside typical range' and AC.assay_id in ('{placeholder}');"

    else:
        assay_query = "select A.chembl_id as assay_chembl_id, A.assay_id, T.chembl_id as target_chembl_id, A.description, A.assay_type, B.site_name, T.target_type, T.pref_name, T.organism,  A.assay_organism, A.assay_strain, A.assay_tissue, A.assay_cell_type, A.confidence_score from assays A left join binding_sites B on A.tid == B.tid left join target_dictionary T on A.tid = T.tid;"
        #cmpd_query = "select md.chembl_id, cs.canonical_smiles, cs.standard_inchi_key from compound_structures cs, molecule_dictionary md where cs.molregno == md.molregno;"
        cmpd_query = "select md.chembl_id as cmpd_chembl_id, cs.canonical_smiles, cp.mw_freebase as mw from compound_structures cs, molecule_dictionary md, compound_properties cp where cs.molregno == md.molregno and cs.molregno == cp.molregno;"
        #activity_query = "select md.chembl_id as cmpd_chembl_id, AY.chembl_id as assay_chembl_id,  AC.standard_value, AC.standard_relation, AC.standard_units, AC.standard_type, AC.data_validity_comment from activities AC, molecule_dictionary md, assays AY where AC.assay_id == AY.assay_id and AC.molregno == md.molregno and AC.data_validity_comment IS NOT 'Outside typical range' and AC.assay_id in ('{placeholder}');"
        
    
    # get All compounds
    cmpd_info = pd.read_sql_query(cmpd_query, conn)
    cmpd_info[['cmpd_chembl_id', 'canonical_smiles']].to_csv(cmpd_sfile, index=False, header=True)
    
    #activity_query = "select md.chembl_id as cmpd_chembl_id, AY.chembl_id as assay_chembl_id,  AC.standard_value as PIC50, AC.standard_relation, AC.standard_units, AC.pchembl_value as pchembl, AC.standard_type from activities AC, molecule_dictionary md, assays AY where AC.assay_id == AY.assay_id and AC.molregno == md.molregno and AC.data_validity_comment IS NOT 'Outside typical range' and AC.assay_id in ({}) and AC.standard_value NOTNULL;".format(','.join([str(x) for x in assay_ids]))
  
    # get all activities
    activity_query = "select md.chembl_id as cmpd_chembl_id, AY.chembl_id as assay_chembl_id, AC.assay_id as assay_id, AC.standard_value as PIC50, AC.standard_relation, AC.standard_units, AC.pchembl_value as pchembl, AC.standard_type from activities AC, molecule_dictionary md, assays AY where AC.assay_id == AY.assay_id and AC.molregno == md.molregno and AC.data_validity_comment IS NOT 'Outside typical range' and AC.standard_value NOTNULL;"
    activities = pd.read_sql_query(activity_query, conn)
    activities_group = activities.groupby('assay_id').count()
    assay_ids = set(activities_group[activities_group.cmpd_chembl_id >= 50].index)
    activities = activities[activities.assay_id.isin(assay_ids)]

    # get assays with at least 50 activities
    if os.path.isfile(assay_file):
        assay_info = pd.read_csv(assay_file, index_col=0, header=0)
    else:
        assay_info = pd.read_sql_query(assay_query, conn)
        assay_info = assay_info[assay_info.assay_id.isin(assay_ids)]
        assay_info = assay_info.groupby('assay_chembl_id').nth(0)
        assay_info.to_csv(assay_sfile, index=True, header=True)

    activities = activities.merge(cmpd_info, how='left', right_on='cmpd_chembl_id', left_on='cmpd_chembl_id')
    del cmpd_info
    activities = activities[['cmpd_chembl_id', 'assay_chembl_id', 'PIC50', 'canonical_smiles', 'standard_relation', 'standard_units', 'pchembl', 'standard_type', 'mw']]
    activities.dropna(axis=0, how='any', subset=['PIC50', 'canonical_smiles'], inplace=True)
    assay_chembl_id_list = list(activities['assay_chembl_id'].unique())

 

    if os.path.isfile(act_sfile):
        df_data = pd.read_csv(act_sfile, index_col=None, header=0, sep='\t')
        existing_assay_chembl_id = list(df_data['assay_chembl_id'].unique())
    else:
        existing_assay_chembl_id = []

    assay_chembl_id_new = [a for a in assay_chembl_id_list if a not in existing_assay_chembl_id]


    output_lock = threading.Lock()
    subLists = partition(assay_chembl_id_new, threads)

    jobs = []

    for i in range(0, threads):

        out_list = subLists[i]
        #thread = threading.Thread(target=ProcessOneByOne, args=(out_list, csv_output_lock))
        thread = multiprocessing.Process(target=query_sql, args=(activities, out_list, preference_first, preference_second, act_sfile, output_lock))
        jobs.append(thread)
    # Start the processes 
    for j in jobs:
        j.start()

    # Ensure all of the processes have finished
    for j in jobs:
        j.join()
    
    print("Finished processing results")
        
    sys.exit(0)

if __name__ == "__main__":
    import os, sys
    import argparse
    import threading
    import multiprocessing
    import sqlite3
    import numpy as np
    import pandas as pd
    main()
