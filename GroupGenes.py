import numpy as np 
import pandas as pd 
import argparse
import os

def generateIndex:
    pass

def main():
    #Parse Args
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--resetIndex', action='store_true',help = 'regenerate the index?')

    CancerTypes = ['BRCA', 'GBM', 'OV', 'LUAD', 'UCEC', 'KIRC',
                    'HNSC', 'LGG', 'THCA', 'LUSC', 'PRAD', 'SKCM'
                    'COAD', 'STAD', 'BLCA', 'LIHC', 'CESC', 'KIRP',
                    'SARC', 'LAML', 'ESCA', 'PAAD', 'PCPG', 'READ',
                    'TGCT', 'THYM', 'THYM', 'KICH', 'ACC', 'MESO',
                    'UVM', 'DLBC', 'UCS', 'CHOL']
    Pipelines = ['Muse', 'Mutect', 'Somatic_Sniper', 'Varscan']
    
    if resetIndex or not os.path.exists()

if __name__ == "__main__":
    main()