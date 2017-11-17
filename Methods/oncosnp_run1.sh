#!/bin/sh

# change path to directory where you have installed OncoSNP or give path to the run_oncosnp.sh script
# batchx$SGE_TASK_ID.txt represents the batch file - split the big batch file into 80 batches each of 5 samples - so we can run OncoSNP in parallel
# the normalcontentlevels are set to 0, 0.3, 0.5 and 0.7

cd PATH_TO_ONCOSNP &&
./run_oncosnp.sh v82/ --batch-file batchx$SGE_TASK_ID.txt --output-dir oncosnp_rez \\
--gcdir b37/ --paramsfile configuration/hyperparameters-affy.dat --levelsfile configuration/levels-affy.dat \\
--trainingstatesfile configuration/trainingStates.dat --tumourstatesfile configuration/tumourStates.dat \\
--intratumour --chr 1 --normalcontentlevels 0:1:0 --allprobes  --fulloutput --plot

