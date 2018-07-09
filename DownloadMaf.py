import numpy as np
import pandas as pd
import argparse

def main():
    CancerTypes = ['BRCA', 'GBM', 'OV', 'LUAD', 'UCEC', 'KIRC',
                    'HNSC', 'LGG', 'THCA', 'LUSC', 'PRAD', 'SKCM'
                    'COAD', 'STAD', 'BLCA', 'LIHC', 'CESC', 'KIRP',
                    'SARC', 'LAML', 'ESCA', 'PAAD', 'PCPG', 'READ',
                    'TGCT', 'THYM', 'THYM', 'KICH', 'ACC', 'MESO',
                    'UVM', 'DLBC', 'UCS', 'CHOL']
    Pipelines = ['Muse', 'Mutect', 'Somatic Sniper', 'Varscan']
    dir = "Mafs/"
    for CancerType in CancerTypes:
        for pipeline in pipelines

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()

    main()
