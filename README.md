Profile-QSAR
============
* Profile-QSAR is project to build massively multitask, two-step machine learning models with unprecedented scope, accuracy, and applicability domain. In step one, a “profile” of conventional single-assay random forest regression models are trained on a very large number of biochemical and cellular pIC50 assays using Morgan 2 substructural fingerprints as compound descriptors. In step two, a panel of partial least squares (PLS) models are built using the profile of pIC50 predictions from those random forest regression models as compound descriptors (hence the name).
* This repo contains (1) scripts for ChEMBL data retrival and profile-QSAR modeling; (2) some example datasets and sample outputs; and (3) one command line tool to look up profile-QSAR pre-calculations. The detailed description of profile-QSAR can be found here: [J Chem Inf Model, 2021](https://doi.org/10.1021/acs.jcim.0c01342); [J Chem Inf Model, 2019](https://doi.org/10.1021/acs.jcim.9b00375); [J Chem Inf Model, 2017](https://doi.org/10.1021/acs.jcim.7b00166); and [J Chem Inf Model, 2011](https://doi.org/10.1021/ci1005004).



The files in this repo:
  * ChEMBL data retrival
    * `./chembl_data_retrival/chembl_sqlite.sh` - a bash script that specify chembl data downloading
    * `./chembl_data_retrival/chembl_sqlite.py` - download assay data, compounds, and meta data from ChEMBL (version 28)

  * profile-QSAR
    * `molproc_pqsar.sh` - a bash script that combines all the model development steps of profile-QSAR together. 
    * `MolProc.py` - a main python script that combines the data processing, deploys the specified python scripts to build models or make predictions, and aggregates intermediate data subsets to master tables.
    * `CommonTools.py` - a core python script that computes morgan fingerprint, splits a whole dataset into training and test subsets, and builds cross-validated partial least squares (PLS) regression models.
    * `ModelBuildingRF.py` - build random forest regression models.
    * `ModelBuildingPLS.py` - build PLS models using a max2 strategy.
    * `ModelPredictions.py` - make predictions on all compounds using random forest or PLS models.
  
  * command line tool
    * `./app/grep_moa.sh` - examples to query profile-QSAR pre-calculations
    * `./app/grepMOA2.py` - look up all the assays hit by a list of compounds greater than a threshold, all the compounds that hit a list of assays greater than a threshold, or all IC50s for a list of compounds and assays.
    * `./app/cid.txt` - sample ChEMBL compounds for query
    * `./app/aid.txt` - sample ChEMBL assays for query
    * `./app/ca_id.txt` - sample ChEMBL compounds and assays for query

  * sample datasets (ChEMBL)
    * `./sample_data/chembl_28_AllAssayData.txt` - sample ChEMBL experimental measurements
    * `./sample_data/chembl_28_AllAssayCmpd.csv` - sample ChEMBL compounds for prediciton
    * `./sample_data/chembl_28_AllAssayInfo.csv` - sample meta-data for ChEMBL assays

  * ChEMBL profile-QSAR models
    * `./chembl_28/RF_Models` - directory for random forest models
    * `./chembl_28/RF_Predictions` - directory for random forest predicitons on all ChEMBL compounds
    * `./chembl_28/PLS_Models` - directory for PLS models
    * `./chembl_28/PLS_Predictions/ZscaledAllAssaysCID.csv` - z-scaled profile-QSAR predicitons (compounds in rows). Z scores were calculated as (y_pred - mean)/stdev. y_pred is the profile-QSAR prediction for on compound. mean and stdev are the average and stdard deviation of predicitons of one particular model for all chembl compounds.   
    * `./chembl_28/PLS_Predictions/ZscaledAllAssaysCID_indices` - a binary file with indexed lines for ZscaledAllAssaysCID.csv
    * `./chembl_28/PLS_Predictions/ZscaledAllAssaysAID.csv` - z-scaled profile-QSAR predicitons (assays in rows)
    * `./chembl_28/PLS_Predictions/ZscaledAllAssaysAID_indices` - a binary file with indexed lines for ZscaledAllAssaysAID.csv
    * `./chembl_28/PLS_Predictions/ZscalingStats.csv` - statistics (e.g., mean and stdev) of profile predicitons by assays
    * `./chembl_28/Summary_pQSAR/SummaryRF.csv` - summary file for random forest models
    * `./chembl_28/Summary_pQSAR/SummaryPLS.csv` - summary file for PLS models
    * `./chembl_28/Summary_pQSAR/SummaryPQSAR.csv` - summary file for profile-QSAR

  * performance test
    * `./tests/pQSAR_tests.ipynb` - a jupyter notebook to show the performance of profile-QSAR



## ChEMBL Data retrival
### Key features
 * Retrieve all the ChEMBL compounds.
 * Retrieve ChEMBL assays with at least 50 measurements. The transformation (e.g., pchembl) is used for most assays. Some unit conversions are performed (e.g., IC50 to pIC50). Assays with standard deviation of >=0.50 will be kept. 
 * Retrieve assays' meta-data, such as target, assay description,	assay_type,	site_name,	target_type,	pref_name, organism.

###  Usage example
 * download the latest SQLite database (chembl_28_sqlite.tar.gz) from this link [chembl database](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
 * unzip the compressed file (tar -xvzf chembl_28_sqlite.tar.gz)
 * edit `chembl_sqlite.sh` to specify `VERSION`, directory of chembl_28.db (`SQLDIRS`), and the destination directory (`DATADIRS`)
 * run the script
```
   $ bash chembl_sqlite.sh

```

## Build profile-QSAR models
### Key features
 * Prepare compounds (e.g., strip ions and standardize smiles) and data for modeling
 * Cluster a dataset and split it into training and test sets. Build single-assay random forest models for all assays
 * Predict the activities of all ChEMBL compounds using those random forest models
 * Build PLS models using random forest predictions for all assays as compound's descriptor. Feature selection (i.e., max2) is applied to select the best models.
 * Predict the activities of all ChEMBL compounds using those PLS models

###  Usage example
* edit `molproc_pqsar.sh` to specify the `VERSION` and destination directory of profile-QSAR models (`DIRS`). Specify the directory of ChEMBL experimental measurements (`PIC50FILE_ORI`),  ChEMBL compounds (`CMPDFILE_ORI`), and meta-data (`AssaysInfo`)
* run `molproc_pqsar.sh`. For example:
```
   $ bash molproc_pqsar.sh

```
## Look up profile-QSAR pre-calculations
### Key features
* Look up all the assays hit by a list of compounds greater than a threshold
* Look up all the compounds that hit a list of assays greater than a threshold
* Look all (p)IC50s for a list of compounds and assays

###  Usage example
`grepMOA2.py` arguments:
```
  -h, --help            show this help message and exit
  -i , --Input          Input file with header line followed by querying compounds' ID and/or assays' AIDs
  -z True/False, --ZpIC50 True/False
                        True (default): Threshold and predictions in Z-scaled values; False: original pIC50 (log molar).
  -t , --Threshold      Threshold to filter out unqualified screening
  -o , --Output         Output file in csv format
  -c, --Compound        Querying by compounds (one compound per row)
  -a, --Assay           Querying by assays (one assay per row)
  -ca, --CA             Querying by compounds(first column) and assays (second column)
```
* look up compounds (CID) with a threshold of 3 stdev above the mean
```
python grepMOA2.py -c -i cid.txt -d ../chembl_28 -t 3.0 -z True -o cid_out.csv
```
	
* look up assays (AID) with a threshold of 3 stdev above the mean
```
python grepMOA2.py -a -i aid.txt -d ../chembl_28 -t 3.0 -z True -o aid_out.csv
```

* look up activities of a list of compounds (CID) and a list of assays (AID)
```
python grepMOA2.py -ca -i ca_id.txt -d ../chembl_28 -z True -o ca_id_out.csv
```

## FAQ
### what is the requirement to run this code?
It requirs a high performance computing cluster to run the code. It has been tested on Novartis High performance Computing cluster using python version 3.8.6 and bash shell. Important packages are in `requirements.txt`.
### Where to edit the code to adapt it to a different cluster?
The python code submits jobs to clusters and gets job ID and status using commands like `qsub`, `subprocess`, and `re`. Users might need to tweak those lines to accommodate the code to their own clusters. 
### Do you fine-tune parameters for larger datasets?
The number of jobs, memoery allocation, R2 cutoff (0.05, 0.20) for max2 have been tuned based on serveral dataset: (1) Novartis (~12K assays and 5.5M compounds) and (2) ChEMBL (a total of 4276 assays and 1.4M compounds) data (refering to [J Chem Inf Model, 2019](https://doi.org/10.1021/acs.jcim.9b00375)). Those same set of paramers also worked well on the example kinase dataset. It might still need to be fine-tuned depending on your own datasets and clusters. 

## Notes and Acknowledgments
@Valery Polyakov, @Li Tian,  @Brian kelly... 

## How do I cite profile-QSAR:
```
@article{Martin2021,
annote = {doi: 10.1021/acs.jcim.0c01342},
author = {Martin, Eric J and Zhu, Xiang-Wei},
doi = {10.1021/acs.jcim.0c01342},
issn = {1549-9596},
journal = {Journal of Chemical Information and Modeling},
month = {apr},
publisher = {American Chemical Society},
title = {{Collaborative Profile-QSAR: A Natural Platform for Building Collaborative Models among Competing Companies}},
url = {https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c01342},
year = {2021}
}

@article{Martin2019,
author = {Martin, Eric J and Polyakov, Valery R and Zhu, Xiang-Wei and Tian, Li and Mukherjee, Prasenjit and Liu, Xin},
doi = {10.1021/acs.jcim.9b00375},
issn = {1549-960X (Electronic)},
journal = {Journal of chemical information and modeling},
language = {eng},
month = {sep},
number = {10},
pages = {4450--4459},
pmid = {31518124},
title = {{All-Assay-Max2 pQSAR: Activity Predictions as Accurate as Four-Concentration IC50s for 8558 Novartis Assays}},
volume = {59},
year = {2019}
}

@article{Martin2017,
author = {Martin, Eric J. and Polyakov, Valery R. and Tian, Li and Perez, Rolando C.},
doi = {10.1021/acs.jcim.7b00166},
issn = {15205142},
journal = {Journal of Chemical Information and Modeling},
number = {8},
pages = {2077--2088},
title = {{Profile-QSAR 2.0: Kinase Virtual Screening Accuracy Comparable to Four-Concentration IC50s for Realistically Novel Compounds}},
volume = {57},
year = {2017}
}

@article{Martin2011,
author = {Martin, Eric and Mukherjee, Prasenjit and Sullivan, David and Jansen, Johanna},
doi = {10.1021/ci1005004},
isbn = {1549-960X (Electronic)\r1549-9596 (Linking)},
issn = {15499596},
journal = {Journal of Chemical Information and Modeling},
number = {8},
pages = {1942--1956},
pmid = {21667971},
title = {{Profile-QSAR: A novel meta-QSAR method that combines activities across the kinase family to accurately predict affinity, selectivity, and cellular activity}},
volume = {51},
year = {2011}
}
```
## Contact Information
For help or issues using profile-QSAR, please submit a GitHub issue.

For personal communication related to this package, please contact Eric Matin (eric.martin@novartis.com) and Xiangwei Zhu (xwzhunc@gmail.com).


