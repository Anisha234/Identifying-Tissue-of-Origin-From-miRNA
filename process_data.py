import pandas as pd
import numpy as np
import glob
import argparse

def assemble_data_frame(df_mRNA,base_dir, columns_to_keep, skiprows=None):
    #this function creates the entire dataframe by concatenating the rows specified in each filename from the info_file
   
    cols_in_total_data=['case_id'] #columns in total_data
    idx=0
    for case_id,id, fname in zip(df_mRNA['cases.0.case_id'], df_mRNA['id'], df_mRNA['file_name']):
        # Load all files for mRNA and filter it

        dir_to_load=base_dir+'\\'+id+'\\'

        file_to_load = dir_to_load+fname
    

        data=pd.read_csv(file_to_load,delimiter='\t',skiprows= skiprows) #data for one person
 
        
        data = data[columns_to_keep]
        if idx==0:
            cols_in_total_data.extend(list(data[columns_to_keep[0]])) #gets the miRNA names. 
            #creates the names for all the columns in an empty dataframe
            total_data=pd.DataFrame(columns=cols_in_total_data) 


        case_data=[case_id]
        #extend the row with the rpms for each sample in for loop
        case_data.extend(list(data[columns_to_keep[1]])) 

        if idx%100 == 0:
            print(idx)
            print(total_data.shape)
        

        total_data.loc[idx]=case_data
        idx=idx+1

    return total_data


def get_meta_data(df,case_id):
    print('In get meta data')
 
  #  df=df.loc[df['cases.0.case_id'].isin(case_id)]
  
    df=df[['cases.0.case_id','cases.0.demographic.gender','cases.0.demographic.days_to_birth','cases.0.project.primary_site','cases.0.samples.0.sample_type']]

    df.replace('male', 1, inplace=True)
    df.replace('female', 0, inplace=True)
    df.replace('Primary Tumor', 0, inplace=True)
    df.replace('Solid Tissue Normal', 1, inplace=True)
    df.replace('Metastatic', 2, inplace=True)
    disease_labels=["Breast", "Uterus", "Ovary", "Prostate", "Testis" ,"Lung", "Kidney","Bladder","Esophagus", "Liver","Pancreas","Pleura","Colorectal", "Skin", "Stomach",
                 "Brain", "Cervix", "Thyroid"]
    df=df.loc[df['cases.0.project.primary_site'].isin(disease_labels)]
  
    df.replace(disease_labels,list(range(len(disease_labels))), #categorical encoding
              inplace=True)
 
    
    df=df.reset_index(drop=True)
    df.rename(columns={'cases.0.case_id':'case_id','cases.0.samples.0.sample_type':'sample_type','cases.0.demographic.gender':'gender','cases.0.demographic.days_to_birth':'age','cases.0.project.primary_site':'disease_type' }, inplace=True)
    df['age'] = -df['age']
    df.replace('False', 0, inplace=True)
    df.replace('not reported',0, inplace=True)
 
    return df

def remove_outliers_zeros(df,  percent_zero_threshold):
    '''
    Only keeps columns where the percent that are zero is < threshold
    '''

    column_cut_off = int(percent_zero_threshold/100*len(df)) 
    b = (df == 0).sum(axis='rows')
    print(b)
    df= df.loc[:, b<=column_cut_off]
    return df

def clean_data(df, percent_zero_threshold):
    print(df.shape)
    df=df.dropna() # drop rows that contain NaNs 
    print(df.shape)
    df=remove_outliers_zeros(df,percent_zero_threshold)
    print(df.shape)
    return df




argParser = argparse.ArgumentParser()
argParser.add_argument('-f', '--info_file',help="file name with file ids")
argParser.add_argument('-d', '--directory', help="directory with miRNA files")
argParser.add_argument("-o", "--output", help="output file name")
args = argParser.parse_args()

df_miRNA =pd.read_csv(args.info_file)
print(df_miRNA.shape)
print(df_miRNA.head())

base_dir= args.directory

total_data_miRNA= assemble_data_frame(df_miRNA,base_dir, columns_to_keep=['miRNA_ID','reads_per_million_miRNA_mapped'])
total_data_miRNA = clean_data(total_data_miRNA, 50)

total_data=total_data_miRNA

total_data.to_csv('pre_metadata_total_data.csv',index=False)
other_data=get_meta_data(df_miRNA, list(total_data['case_id']))
print(other_data.shape)
other_data.to_csv('other_data.csv',index=False)

print('Otherdata', other_data.head())
print('Totaldata', total_data.head())
# Now add gender and age values

total_data['age'] = other_data['age']
total_data['gender'] = other_data['gender']
total_data['disease_type'] = other_data['disease_type']
total_data['sample_type'] = other_data['sample_type']
print(total_data.head())
total_data['age'].fillna((total_data['age'].mean()), inplace=True)
total_data['gender'].fillna(False,inplace=True)
total_data['age'] = (total_data['age'] - total_data['age'].mean())/total_data['age'].std()
print(total_data.head())
print(total_data.shape)
total_data.to_csv(args.output,index=False)