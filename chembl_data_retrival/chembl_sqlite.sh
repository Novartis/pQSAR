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

VERSION=chembl_28
# directory of chembl_28.db
SQLDIRS=chembl_28.db
DATADIRS=./sample_data
threads=40

python chembl_sqlite.py ${SQLDIRS} ${DATADIRS} ${VERSION} False ${threads}

