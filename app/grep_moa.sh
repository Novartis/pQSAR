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

#look up compounds (CID) with a threshold of 3 stdev above the mean
python grepMOA2.py -c -i cid.txt -d ../chembl_28 -t 3.0 -z True -o cid_out.csv

#look up assays (AID) with a threshold of 3 stdev above the mean
python grepMOA2.py -a -i aid.txt -d ../chembl_28 -t 3.0 -z True -o aid_out.csv

#look up activities of a list of compounds (CID) and a list of assays (AID)
python grepMOA2.py -ca -i ca_id.txt -d ../chembl_28 -z True -o ca_id_out.csv
