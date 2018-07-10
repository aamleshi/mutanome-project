import numpy as np
import pandas as pd
import argparse
import os
from tqdm import tqdm
import requests
import json
import re
import gzip

def check_gdc_status(filepath):
    '''makes the request to gdc server to download a specific cancer type'''
    print(filepath)
    status_endpt = "https://api.gdc.cancer.gov/status"
    response = requests.get(status_endpt)
    print(response.content)

def downloadMaf(filepath, cancerType, pipeline):
    
    if pipeline == "Muse":
        workflow = "MuSE Variant Aggregation and Masking"
    elif pipeline == 'Mutect':
        workflow = "MuTect2 Variant Aggregation and Masking"
    elif pipeline == 'Somatic_Sniper':
        workflow = "SomaticSniper Variant Aggregation and Masking"
    elif pipeline == 'Varscan':
        workflow = "VarScan2 Variant Aggregation and Masking"


    
    files_endpt = "https://api.gdc.cancer.gov/files"

    
    filters = {
        "op": "and",
        "content":[
            {
            "op": "in",
            "content":{
                "field": "data_format",
                "value": ["MAF"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.cases.project.project_id",
                "value": ["TCGA-" +cancerType]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "access",
                "value": ["open"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.analysis.workflow_type",
                "value": [workflow]
                }
            }
        ]
    }

    print("TCGA-" +cancerType)

    # Here a GET is used, so the filter parameters should be passed as a JSON string.

    params = {
        "filters": json.dumps(filters),
        "fields": "file_id",
        "format": "JSON",
        "size": "100"
        }

    response = requests.get(files_endpt, params = params)

    print(response.content)
    
    file_uuid_list = []

    # This step populates the download list with the file_ids from the previous query
    for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
        file_uuid_list.append(file_entry["file_id"])

    data_endpt = "https://api.gdc.cancer.gov/data"

    params = {"ids": file_uuid_list}

    response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

    response_head_cd = response.headers["Content-Disposition"]

    #file_name = re.findall("filename=(.+)", response_head_cd)[0]

    with open(filepath+".gz", "wb") as output_file:
        output_file.write(response.content)
    zipedFile = gzip.open(filepath+".gz", 'rb')
    f = open(filepath, 'wb')
    f.write(zipedFile.read())
    zipedFile.close()
    f.close()
    os.remove(filepath+".gz"l)

def main():
    #Parse Args
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--reset', action='store_true',help = 'redownload existing files?')
    #Cancer type and pipeline data, maybe allow user to only download a subset?
    CancerTypes = ['BRCA', 'GBM', 'OV', 'LUAD', 'UCEC', 'KIRC',
                    'HNSC', 'LGG', 'THCA', 'LUSC', 'PRAD', 'SKCM'
                    'COAD', 'STAD', 'BLCA', 'LIHC', 'CESC', 'KIRP',
                    'SARC', 'LAML', 'ESCA', 'PAAD', 'PCPG', 'READ',
                    'TGCT', 'THYM', 'THYM', 'KICH', 'ACC', 'MESO',
                    'UVM', 'DLBC', 'UCS', 'CHOL']
    print(len(CancerTypes))
    Pipelines = ['Muse', 'Mutect', 'Somatic_Sniper', 'Varscan']

    top_dir = "MafArchive/"
    if not os.path.exists(top_dir):
        os.mkdir(top_dir)
    for cancerType in tqdm(CancerTypes):
        cancer_dir = top_dir + cancerType + '/'
        if not os.path.exists(cancer_dir):
            os.mkdir(cancer_dir)
        for pipeline in Pipelines:
            filepath = cancer_dir+pipeline+'.maf'
            if not os.path.exists(filepath):
                downloadMaf(filepath, cancerType, pipeline)

if __name__ == "__main__":
    main()
