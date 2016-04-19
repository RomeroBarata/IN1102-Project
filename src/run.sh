#!/bin/bash

# This script runs 'mvfcmv' using the distance matrices required. The
# user must provide the following arguments:
#   - the path containing the matrices to be used;

DATA_PATH=$1
# q=1*, k=10*, m=2*, epsilon=1e-10*, 4 instances, R data into
# 'mvfcmv-data.R'
# * - fixed values for this project.
PRG_ARGS="-q 1 -k 10 -m 2 -e 1e-10 -i 4 --Rfile mvfcmv-data.R"

./mvfcmv $PRG_ARGS 2000 3 $DATA_PATH/mfeat-fac-diss $DATA_PATH/mfeat-fou-diss $DATA_PATH/mfeat-kar-diss
