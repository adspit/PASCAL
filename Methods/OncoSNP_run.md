For OncoSNP, we first downloaded and installed the algorithm by following the instructions given at: https://sites.google.com/site/oncosnp/.
Once the installation was done, we ran the following command:
```
./run_oncosnp.sh  --batch-file batch_file.txt --output-dir output_oncosnp \
  --gcdir b37 --hgtables configuration/hgTables_b37.txt --paramsfile configuration/hyperparameters-affy.dat --levelsfile configuration/levels-affy.dat \
  --trainingstatesfile configuration/trainingStates.dat --tumourstatesfile configuration/tumourStates.dat --isdiploid --chr 1 --fulloutput --allprobes

´´´
