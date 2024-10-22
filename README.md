#Description of folder
+ dataSet
  + density (fold): Dataset selected using density.
  + pc (fold): Dataset selected using packing coefficient.
  + D_chon.csv: Dataset for CHON molecules.
  + D_no2.csv:Dataset for nitro-containing molecules.
  + data_run.py: Data processing code, including test set construction, farthest point sampling of molecular fingerprints, data cleaning, etc.
  + testSet.csv: Independent test set
+ model 
  + generalizability_test:Effect of the difference between test set and training set on the generalization ability of the models.
  + model_density
  + model_D_chon
  + model_D_no2
  + model_D_pc  
  
  
All the models were trained with [chemProp](https://github.com/chemprop/chemprop), detailed parameters see the args.json file of each model.

