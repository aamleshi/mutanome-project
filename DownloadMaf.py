import numpy as np
import pandas as pd
import argparse
import os
from tqdm import tqdm
import requests
import json


def downloadMaf(filepath, cancerType, pipeline):
    '''makes the request to gdc server to download a specific cancer type'''
    print(filepath)
    status_endpt = "https://api.gdc.cancer.gov/status"
    response = requests.get(status_endpt)
    print(response.content)


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
    Pipelines = ['Muse', 'Mutect', 'Somatic_Sniper', 'Varscan']
    ####DEBUG LIMITER####
    CancerTypes = CancerTypes[0:1]
    ####DEBUG LIMITER####
    top_dir = "MafArchive/"
    if not os.path.exists(top_dir):
        os.mkdir(top_dir)
    for cancerType in tqdm(CancerTypes):
        cancer_dir = top_dir + cancerType + '/'
        if not os.path.exists(cancer_dir):
            os.mkdir(cancer_dir)
        for pipeline in Pipelines:
            filepath = cancer_dir+pipeline
            if not os.path.exists(filepath):
                downloadMaf(filepath, cancerType, pipeline)

if __name__ == "__main__":
    main()
