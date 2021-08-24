#!/bin/sh
#"""
#Copyright 2021 Novartis Institutes for BioMedical Research Inc.
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#"""


########################
# CODEDIRS: the directory with all the necessary python code:
# MDLDIRS: the directory to save all the pQSAR files
# PIC50FILE: tab separated file with a header
	# compoundID	AssayID	PIC50	smiles
# CMPDFILE: csv or txt file with a header
	# compoundID,smiles
########################

VERSION=chembl_28

CODEDIRS=$PWD
DIRS=./${VERSION}

PIC50FILE_ORI=./sample_data/${VERSION}_AllAssayData.txt
CMPDFILE_ORI=./sample_data/${VERSION}_AllAssayCmpd.csv
AssaysInfo=./sample_data/${VERSION}_AllAssayInfo.csv

PIC50FILE=${DIRS}/${VERSION}_AllAssayData.txt
CMPDFILE=${DIRS}/${VERSION}_AllAssayCmpd.csv

mkdir $DIRS
MDLDIRS=`realpath $DIRS`

#######################

time python ${CODEDIRS}/MolProc.py -j DataPrep -i ${PIC50FILE_ORI} -f ${CMPDFILE_ORI} -f2 ${PIC50FILE} -o ${CMPDFILE} -p ${CODEDIRS}/DataPreparation.py

# RF models building
time python ${CODEDIRS}/MolProc.py -j RF4pQSAR -o ${MDLDIRS} -i ${PIC50FILE} -p ${CODEDIRS}/ModelBuildingRF.py 

# RF prediction
time python ${CODEDIRS}/MolProc.py -j RFpredict -o ${MDLDIRS} -d ${MDLDIRS}/RF_Models -i ${CMPDFILE} -p ${CODEDIRS}/ModelPredictions.py 

# PLS models building
time python ${CODEDIRS}/MolProc.py -j PLS4pQSAR -o ${MDLDIRS} -i ${MDLDIRS}/RF_Models/PredvsExp.csv -d ${MDLDIRS}/RF_Predictions -p ${CODEDIRS}/ModelBuildingPLS.py

# PLS prediction
time python ${CODEDIRS}/MolProc.py -j PLSpredict -o ${MDLDIRS} -d ${MDLDIRS}/PLS_Models -i ${MDLDIRS}/RF_Predictions/MasterTable.csv -cpd ${MDLDIRS}/Summary_pQSAR/CIDReference.csv -p ${CODEDIRS}/ModelPredictions.py 

#Create summary table
time python ${CODEDIRS}/MolProc.py -j pQSARsummary -o ${MDLDIRS}/Summary_pQSAR -i ${MDLDIRS}/Summary_pQSAR/SummaryRF.csv -f ${MDLDIRS}/Summary_pQSAR/SummaryPLS.csv -f2 ${MDLDIRS}/PLS_Predictions/ZscalingStats.csv -f3 ${AssaysInfo} -f4 ${PIC50FILE}
