import numpy as np 
import pandas as pd 
import argparse
import os
from tqdm import tqdm
import pickle
import requests

def createGeneIndex(mafDir, cancerTypes, pipelines, geneIndexFile):
    #Creates a dataframe that has all unique genes
    #importaint column is the protein files, which will be populated later
    geneFrame = pd.DataFrame(columns = ['Entrez_Gene_Id', 'uniProt' 'pdb_files']+pipelines)
    for cancer in tqdm(cancerTypes):
        cancerDir = mafDir+ cancerType + '/'
        for pipeline in pipelines:
            filepath = cancerDir+pipeline+'.maf'
            with open(filepath, r) as mafFile:
                mafDF = pd.read_table(file, skiprows= 5)
                geneIds = mafDF['Entrez_Gene_Id'].unique()
                for geneId in geneIds:
                    try:
                        geneFrame.loc[geneId, pipeline] = 1
                    except:
                        geneFrame = geneFrame.append({'Entrez_Gene_Id': 1, pipeline:1}, ignore_index=True)

    geneFrame = geneFrame.set_index('Entrez_Gene_Id')
    #now pickle the geneFrame
    geneFrame.to_pickle(geneIndexFile)

def getUniprot(Entrez, split = True):
    #Takes in an entrez id, and receives the top uniprot id
    """Returns the first Uniprot named structure for a entrez gene id"""
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from':'P_ENTREZGENEID',
    'to':'ID',
    'format':'tab',
    'query':Entrez
    }

    data = urllib.parse.urlencode(params)
    request = urllib.request.Request(url, data.encode("utf-8"))
    contact = "aa@uchicago.edu" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(request)
    page = response.read(200000)
    page = page.decode("utf-8")
    if split:
        try:
            return page.split('\n')[1].split('\t')[1]
        except:
            return page
    else:
        return page

def mapUniProt(geneIndexFile):
    #fills in the uniProt column of the geneIndex dataframe by calling the uniprot mapping server
    tqdm.pandas(desc="mapping entrez to uniProt")
    indexFrame = pd.read_pickle(geneIndexFile)
    indexFrame['uniProt'] = indexFrame['Entrez_Gene_Id'].progress_apply(lambda entrez: getUniProt(entrez))

def associateGenePDB(geneIndexFile):
    #associates the entrezGeneIds with a sequence of pdb files
    #use a greedy algorithm to determine which sequences to use
    pass


def main():
    #Parse Args
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--resetIndex', action='store_true',help = 'regenerate the index?')
    arg_parser.add_argument('--resetPdb', action='store_true',help = 'ping servers to re-associate geneIds with pdb files')

    cancerTypes = ['BRCA', 'GBM', 'OV', 'LUAD', 'UCEC', 'KIRC',
                    'HNSC', 'LGG', 'THCA', 'LUSC', 'PRAD', 'SKCM'
                    'COAD', 'STAD', 'BLCA', 'LIHC', 'CESC', 'KIRP',
                    'SARC', 'LAML', 'ESCA', 'PAAD', 'PCPG', 'READ',
                    'TGCT', 'THYM', 'THYM', 'KICH', 'ACC', 'MESO',
                    'UVM', 'DLBC', 'UCS', 'CHOL']
    pipelines = ['Muse', 'Mutect', 'Somatic_Sniper', 'Varscan']
    mafDir = "MafArchive/"
    geneIndexFile = "geneIndex.pkl"
    if resetIndex or not os.path.exists(geneIndexFile):
        createGeneIndex(mafDir, cancerTypes, pipelines, geneIndexFile)
        associateGenePDB(geneIndexFile)
    if resetPdb:
        associateGenePDB(geneIndexFile)
    
    

if __name__ == "__main__":
    main()
