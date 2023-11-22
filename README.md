# miRNA

This repository has four files in conjunction with our research paper. The files are gdc_query.py, process_data.py, SRA_preprocess.ipynb and miRNA_model_training_eval_new.ipynb. 
These files do the following:
1. gdc_query.py: This file queries data from the GDC data portal to obtain samples from the TCGA repository.
2. process_data.py: This file is responsible for data processing of the TCGA file. It produces a single dataframe that can be used for machine learning experiments.
3. SRA_preprocess.ipynb: This notebook preprocesses the data downloaded from the MITed portal (corresponding to the SRA dataset). Specifically, this notebook identifies common micro RNAs between TCGA and SRA.
4. miRNA_model_training_eval_new.ipynb: This notebook contains all the classifiers and the results of our machine learning experiments. 
