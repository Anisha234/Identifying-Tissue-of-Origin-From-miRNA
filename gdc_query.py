import requests
import json
import re
import pandas as pd
import argparse

fields = [
    "file_name",
    "cases.case_id",
    "cases.samples.sample_type",
    'cases.primary_site',
    'cases.project.primary_site',
    "cases.demographic.gender",
    "cases.demographic.days_to_birth",
    #"cases.exposures.years_smoked",
    #"cases.exposures.alcohol_history",
    ]

fields = ",".join(fields)

files_endpt = "https://api.gdc.cancer.gov/files"

def get_filter(primary_site, tissue_type, mode):
# This set of filters is nested under an 'and' operator.
    if mode =='mRNA':
        data_category = "transcriptome profiling"
        data_type="Gene Expression Quantification"
        data_format="TSV"
    if mode =='miRNA':
        data_category = "transcriptome profiling"
        data_type="miRNA Expression Quantification"
        data_format="TXT"
    if mode=='Methylation':
        data_category='DNA Methylation'
        data_type="Methylation Beta Value"
        data_format="TXT"
 

    filters = {
        "op": "and",
        "content":[
            {
            "op": "in",
            "content":{
                "field": "cases.project.primary_site",
                "value": primary_site
                }
            },
            {
            "op": "in",
            "content":{
                "field": "cases.samples.sample_type",
                "value": tissue_type
              #  "value": ["Blood Derived Normal", "Solid Tissue Normal"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_category",
                "value": [data_category]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_type",
                "value": [data_type]
                }
            },
           {
            "op": "in",
            "content":{
                "field": "files.data_format",
                "value": [data_format]
                }
            }
        ]
    }
    return filters


argParser = argparse.ArgumentParser()
argParser.add_argument('-n', '--cancer_list', nargs='+', default=[])
argParser.add_argument('-t', '--tissue_type', nargs='+', default=[])
argParser.add_argument("-m", "--data", help="data (mRNA, miRNA, Methylation)")
args = argParser.parse_args()


mode = args.data
cancer_site=args.cancer_list
print(cancer_site)
filters=get_filter(cancer_site, args.tissue_type, mode)

# A POST is used, so the filter parameters can be passed directly as a Dict object.
params = {
    "filters": filters,
    "fields": fields,
    "format": "CSV",
    "size": "10000"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)
res=response.content.decode("utf-8")
print(type(res))


with open("my_file_met.csv", "wb") as binary_file:
   
    # Write bytes to file
    binary_file.write(response.content)


  
# assign dataset
csvData = pd.read_csv("my_file_met.csv")
                                         
# displaying unsorted data frame
print("\nBefore sorting:")
print(csvData)
  
# sort data frame
csvData.sort_values(["cases.0.case_id"], 
                    axis=0,
                    ascending=[False], 
                    inplace=True)
  
# displaying sorted data frame
print("\nAfter sorting:")
print(csvData)



# Remove rows where a key occurs only once
#fname = 'df2_'+ mode +'.csv'
#csvData.to_csv(fname,index=False)

def download_data(file_ids):
# Now get the file ids

    data_endpt = "https://api.gdc.cancer.gov/data"

    params = {"ids": file_ids}

    response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})
 
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]
    print(file_name)





    with open(file_name, "wb") as output_file:
        output_file.write(response.content)
    return

file_ids = list(csvData['id'])
download_data(file_ids)

