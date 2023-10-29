import pandas as pd
import numpy as np
import glob
import argparse

def assemble_data_frame(mode, df_mRNA, c_mRNA,base_dir, columns_to_keep, extension='tsv',skiprows=None):
    
    total_data=pd.DataFrame()
    genes=['case_id']
    print('Number of Files', len(c_mRNA))
    idx=0
    for case_id,id, fname in zip(df_mRNA['cases.0.case_id'], df_mRNA['id'], df_mRNA['file_name']):
        # Load all files for mRNA and filter it

        dir_to_load=base_dir+'\\'+id+'\\'
       
        file_to_load = dir_to_load+fname
       
        if mode =='Methylation':
           data=pd.read_csv(file_to_load,delimiter='\t',skiprows= skiprows,header=None, names=['gene_name', 'Methylation'])
        else:
           data=pd.read_csv(file_to_load,delimiter='\t',skiprows= skiprows)
 
   
        genes.extend(list(data[columns_to_keep[0]]))
        data = data[columns_to_keep]


    
        if idx==0:
            total_data=pd.DataFrame(columns=genes)
            print('idx0', total_data)

     
        case_data=[case_id]
        case_data.extend(list(data[columns_to_keep[1]]))
     
        if idx%100 == 0:
            print(idx)
            print(total_data.shape)
        if mode=='Methylation':
            newFrame=pd.DataFrame(columns=genes)
            newFrame.loc[0]=case_data
            total_data=total_data.merge(newFrame,how='outer')
            print('TotalData', total_data.shape)
        else:
            total_data.loc[idx]=case_data
        idx=idx+1

    return total_data

def assemble_data_frame_old(mode, df_mRNA, c_mRNA,base_dir, columns_to_keep, extension='tsv',skiprows=None):
    
    total_data=pd.DataFrame()
    genes=['case_id']
    print('Number of Files', len(c_mRNA))
    for idx in range(len(c_mRNA)):
        # Load all files for mRNA and filter it
        row= df_mRNA.loc[df_mRNA['cases.0.case_id'] == c_mRNA[idx]]
        # Take first of multiple matches
        fid=list(row['id'])
        case_id=c_mRNA[idx]
        dir_to_load=base_dir+'\\'+str(fid[0])+'\\'
     
        if mode =='Methylation':
           data=pd.read_csv(file_to_load[0],delimiter='\t',skiprows= skiprows,header=None, names=['gene_name', 'Methylation'])
        else:
           data=pd.read_csv(file_to_load[0],delimiter='\t',skiprows= skiprows)
 
      
        genes.extend(list(data[columns_to_keep[0]]))
        data = data[columns_to_keep]
   
        if idx==0:
            total_data=pd.DataFrame(columns=genes)
            print('idx0', total_data)

      
        case_data=[case_id]
        case_data.extend(list(data[columns_to_keep[1]]))
      
        if idx%100 == 0:
            print(idx)
            print(total_data.shape)
        if mode=='Methylation':
            newFrame=pd.DataFrame(columns=genes)
            newFrame.loc[0]=case_data
            total_data=total_data.merge(newFrame,how='outer')
            print('TotalData', total_data.shape)
        else:
            total_data.loc[idx]=case_data
   
    return total_data

def remove_outliers(dataframe,  threshold_low, threshold_hi):
    '''
    Function to remove outliers where less than a threshold and higher than another threshold from a dataframe
    Returns the updated dataframe with the outliers removed.
    '''
    df1 = dataframe.sum(axis=1)
    dataframe['sum'] = df1

    low_valid = dataframe['sum'] > threshold_low
    hi_valid = dataframe['sum'] < threshold_hi
    column_valid = [a and b for a,b in zip(low_valid,hi_valid)] #take all columns satisfiying low valid and high valid

    dataframe = dataframe.drop(columns=['sum'])
    truncated_data = dataframe[column_valid]
    removeAmount = list(dataframe.shape)[0] - list(truncated_data.shape)[0]

    return truncated_data

def get_meta_data(df,case_id):
    print('In get meta data')
 

    df=df[['cases.0.case_id','cases.0.demographic.gender','cases.0.demographic.days_to_birth','cases.0.project.primary_site','cases.0.samples.0.sample_type']]

    df.replace('male', 1, inplace=True)
    df.replace('female', 0, inplace=True)
    df.replace('Primary Tumor', 0, inplace=True)
    df.replace('Solid Tissue Normal', 1, inplace=True)
    df.replace('Metastatic', 2, inplace=True)
    disease_labels=["Breast", "Uterus", "Ovary", "Prostate", "Testis" ,"Lung", "Kidney","Bladder","Esophagus", "Liver","Pancreas","Pleura","Colorectal", "Skin", "Stomach",
                 "Brain", "Cervix", "Thyroid"]
    df=df.loc[df['cases.0.project.primary_site'].isin(disease_labels)]
  
    df.replace(disease_labels,list(range(len(disease_labels))),
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
    df=df.dropna()
    print(df.shape)
    df=remove_outliers_zeros(df,percent_zero_threshold)
    print(df.shape)
    return df

def normalize_data(df):
    print('df', df)
    data_part = df.loc[:, df.columns != 'case_id']
    normalized_df=(data_part-data_part.mean(numeric_only=True))/data_part.std(numeric_only=True)

    normalized_df.insert(0, 'case_id', df['case_id'])
    normalized_df['case_id'] = df['case_id']
    print(normalized_df.head())
    return normalized_df



argParser = argparse.ArgumentParser()
argParser.add_argument('-f', '--info_file',help="file name with file ids")
argParser.add_argument('-d', '--directory', help="directory with miRNA files")
argParser.add_argument("-o", "--output", help="output file name")
args = argParser.parse_args()

df_miRNA =pd.read_csv(args.info_file)
print(df_miRNA.shape)
print(df_miRNA.head())

c_miRNA = list(df_miRNA['cases.0.case_id'])

print(len(c_miRNA))

print(df_miRNA.head())
base_dir= args.directory




total_data_miRNA= assemble_data_frame('miRNA', df_miRNA, c_miRNA,base_dir, columns_to_keep=['miRNA_ID','reads_per_million_miRNA_mapped'], extension='txt')
total_data_miRNA = clean_data(total_data_miRNA, 50)
total_data_miRNA=normalize_data(total_data_miRNA)

# Now concatenate and write the total data [case_id, <all data>, gender, age, disease_label]

total_data=pd.DataFrame()
total_data=total_data_miRNA



total_data.to_csv('pre_merge_total_data.csv',index=False)
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