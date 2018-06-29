import re
import gzip
import glob
import os.path
import urllib
import urllib2
import httplib
import pprint
import string
import math
import random
import copy
from classes.protein_info import PDB
from Bio import SwissProt
from Bio import SeqIO
from Bio import pairwise2
from Bio.PDB import *
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1


import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import pymol

def read_pdb(pdbID, source):

  # Check and download biological assemblies
  assemfound = 0
  if source == "ASSEMBLY" and len(pdbID) < 5: 

      pdbpath ='pdb/assemblies/' 
      pdbfile = pdbID.lower() + '.pdb'

      # Check if biological assembly files already exist
      if glob.glob(pdbpath + pdbfile + '*'):
          assemfound = 1       
      else:
          # Looking for assemblies
          url = 'http://www.rcsb.org/pdb/rest/getEntityInfo?structureId=' + pdbID
          req = urllib2.Request(url)
          f = urllib2.urlopen(req)
          result = f.read()
          regassem = re.search("bioAssemblies=\"(\d)\"", result)

          if hasattr(regassem, 'group'):
              nassem = int(regassem.group(1))

              # Download assemblies
              for a in range(1,nassem+1):
        
                 # Check if file exists locally
                 localfile = pdbpath + pdbfile + str(a)        
                                    
                 # Check if resource exists
                 try:
                     check = urllib2.urlopen("http://www.rcsb.org/pdb/files/" + pdbfile + str(a) + '.gz') 
                     status = check.code
                 except urllib2.HTTPError, e:
                     status = e.code                

                 # Download if resource exists
                 if status == 200:
                    assemfound = 1
             
                    # Download biological assembly
                    print "   Downloading " + pdbID + " biological assembly " + str(a) + " ..."
                    url = 'http://www.rcsb.org/pdb/files/' + pdbfile + str(a) + '.gz'
                    urllib.urlretrieve(url, localfile + ".gz")

                    # Unzip biological assembly
                    gz = gzip.open(localfile + ".gz", 'rb')
                    with open(localfile, 'wb') as out: 
                       out.writelines(gz) 
                    gz.close() 
                    os.remove(localfile + ".gz") 

  # Check and donwload if source is PDB or no biological assemblies have been found
  elif source == "PDB" or assemfound == 0:

      # Local path for PDB
      pdbfile = pdbID.lower() + ".pdb"
      pdbpath = "pdb/pdb/"
  
      # Download file if it is not exist
      if not os.path.exists(pdbpath + pdbfile):
          print "   Downloading " + pdbID + " file..."
          if len(pdbID) < 5:
             url = "http://files.rcsb.org/download/" + pdbfile
          else:
             url = "https://swissmodel.expasy.org/repository/uniprot/" + pdbfile + "?provider=swissmodel"   
          f = urllib2.urlopen(url)
          data = f.read()
          with open(pdbpath + pdbfile, "wb") as code:
              code.write(data)  

  # Read structure from file 
  pathfiles = glob.glob(pdbpath + pdbfile + '*')  
  if pathfiles and not os.stat(pathfiles[0]).st_size == 0:    
     structures = []
     parser = PDBParser(QUIET=True)
     for fname in pathfiles:
        pdbname = fname.split("/")[-1]
        pdbname = pdbname.split(".pdb")
        pdbname = [x for x in pdbname if x]
        pdbname = '.'.join(pdbname)
        structures += [parser.get_structure(pdbname, fname)]
     return structures
  else:
     return None

def pdb_distances(pdb, pdb_seq_pos, outpath, pdbsource):

   pymol.finish_launching()

   # Read and format inputs
   chain = pdb.chain
   var_pdb_pos = pdb.var_pdb_pos
   var_up_pos = sorted(pdb.var_up_pos)
   var_freqs = pdb.freq

   # Recover file name
   filename = pdb.name.split('.')
   struct = filename[0]
   if len(filename)>1:
       filename = '.pdb'.join(filename)
   else:
       filename = struct + '.pdb'

   # Find path for structure
   if pdbsource == "PDB":
       structPath = "pdb/pdb/" + filename
   elif pdbsource == "MODBASE":
       structPath = "modbase/" + filename
       chain = ""
   elif pdbsource == "ASSEMBLY":
       structPath = "pdb/assemblies/" + filename
       if not os.path.exists(structPath):
           structPath = "pdb/pdb/" + filename

   # Load structure to PyMol
   pymol.cmd.load(structPath, struct)

   # Calculate distances
   numvar = len(var_pdb_pos)
   distances  = [[0 for x in range(numvar)] for y in range(numvar)]    
   waps = [[0 for x in range(numvar)] for y in range(numvar)]

   #random_vars = []
   #random_rep = 1000
   soft_cutoff = 6
   #random_wap = [0] * random_rep
   #for x in range(0,random_rep):
   #   random_vars.append(random.sample(pdb_seq_pos, numvar))   

   print "Calculating distances for " + pdb.name + " ..."

   for i in range(0,numvar):     
      # Format variant 1
      pos_init = "/" + struct + "//" + chain + "/{0}".format(var_pdb_pos[i])      
      for j in range(i+1, numvar):
          # Format variant 2
          pos_end = "/" + struct + "//" + chain + "/{0}".format(var_pdb_pos[j])          

          # Calculate distance   
          distance = pymol.cmd.distance('tmp',pos_init,pos_end,mode=4)
          distances[i][j] = distance
          distances[j][i] = distance
          
          # Calculate wap
          norm_freq_i = float(var_freqs[i]**3) / float(2**3 + var_freqs[i]**3)
          norm_freq_j = float(var_freqs[j]**3) / float(2**3 + var_freqs[j]**3)
          waps[i][j] = norm_freq_i * norm_freq_j * math.exp(-float(distance**2) / float((2*(soft_cutoff**2))))  

   # Save distances in a file
   distance_file = outpath + 'distances/' + pdb.name + '_distances.txt'
   fdispdb = open(distance_file, 'w')
   fdispdb.write("Pos:\t" + '\t'.join(str(x) for x in var_up_pos) + '\n')
   posd = 0
   for d in distances:
       fdispdb.write(str(var_up_pos[posd]) + "\t" + '\t'.join(str(x) for x in d) + '\n')
       posd += 1 
   fdispdb.close()

   return (distances, waps)


def check_pdb_coverage(pdbobj, UPseq, structure, outpath, pdbsource):
 
   pdb = copy.copy(pdbobj)

   # Read information from pdb object
   pdbID = pdb.name
   chain = pdb.chain

   # Return if chain is not specified
   if chain == '@':
      return None

   # Check included chains
   chains = chain.split("/")
   print "CHAINS: " 
   incchains = [ch.id for ch in structure[0] if ch.id in chains]
   print incchains
   if not incchains:
      return None
   else:
      chain = incchains[0]

   PDBseq = ""
   pdbPositions = []
   for res in structure[0][chain]:
       if is_aa(res):
           pdbPositions.append(res.id[1])
           PDBseq = PDBseq + seq1(res.resname)        

   # Convert sequence, variants positions and variants AAs to UP-PDB format  
   UP_PDB_pos = [x - pdb.start for x in pdb.var_up_pos]
   UP_PDBseq = UPseq[pdb.start-1:pdb.end]
   UP_PDBseq = string.replace(UP_PDBseq,'U','X')

   # Align UP-PDB and PDB sequences
   alignment = pairwise2.align.globalds(UP_PDBseq, PDBseq, matlist.blosum62, -20, -0.5)

   print alignment

   # Index all AAs in alignment
   seqUP_aa_index = [i for i, x in enumerate(alignment[0][0]) if x != '-']
   seqPDB_aa_index = [i for i, x in enumerate(alignment[0][1]) if x != '-']
   
   # Retrieve variants positions and AAs in alignment
   aln_var_pos = [seqUP_aa_index[i] for i in UP_PDB_pos if alignment[0][1][seqUP_aa_index[i]] != '-']
   alnUP_var = [alignment[0][0][i] for i in aln_var_pos]
   alnPDB_var = [alignment[0][1][i] for i in aln_var_pos]
   if '-' in alnPDB_var:
       print "In here"
       return None

   # Convert variant positions and AAs to PDB sequence  
   seqPDB_pos = [i for i, x in enumerate(seqPDB_aa_index) if x in aln_var_pos]
   seqUP_pos = [i+pdb.start for i, x in enumerate(seqUP_aa_index) if x in aln_var_pos]
   seqPDB_var_pos = [pdbPositions[i] for i in seqPDB_pos]
   seqPDB_var = [PDBseq[i] for i in seqPDB_pos]
 
   # Sorted variants according to the variant position
   sorted_idx = [i[0] for i in sorted(enumerate(pdb.var_up_pos), key=lambda x:x[1])]
   sorted_variants = []
   sorted_positions = []
   sorted_projects = []
   sorted_freqs = []
   for i in sorted_idx:
       if pdb.var_up_pos[i] in seqUP_pos:  
          sorted_variants += [pdb.variants[i]]
          sorted_positions += [pdb.var_up_pos[i]]
          sorted_projects += [pdb.project[i]]
          sorted_freqs += [pdb.freq[i]]

   # Update PDB object with updated variant positions
   pdb.var_pdb_pos = seqPDB_var_pos
   pdb.var_up_pos = sorted_positions
   pdb.variants = sorted_variants
   pdb.project = sorted_projects
   pdb.freq = sorted_freqs
   pdb.num_variants = len(pdb.var_pdb_pos)
   pdb.name = structure.id
   pdb.chain = chain

   # Calculate distances to the updated PDB
   distances, wap = pdb_distances(pdb, pdbPositions, outpath, pdbsource)
   pdb.distances = distances
   pdb.wap = wap

   return pdb

def create_clusters(pdb, cutoff):

   # Retrieve list of variants 
   var_pdb_pos = pdb.var_pdb_pos
   var_up_pos = sorted(pdb.var_up_pos)
   distances = pdb.distances
   waps = pdb.wap

   # Number of variants
   numvar = len(distances)

   # Create cluster 
   clusters_up = []
   clusters_pdb = []   
   clusters_wap = []
   for i in range(0,numvar):
        for j in range(i+1,numvar):
            added = 0
            if distances[i][j] <= cutoff:

                 if not clusters_up:
                    clusters_up.append([var_up_pos[i], var_up_pos[j]]) 
                    clusters_pdb.append([var_pdb_pos[i], var_pdb_pos[j]])             
                    clusters_wap.append(waps[i][j])
                 else:
                    for cl in range(0, len(clusters_up)):
                        c = clusters_up[cl] 
                        if var_up_pos[i] in c and var_up_pos[j] in c:
                            added = 1
                            continue
                        elif var_up_pos[i] in c:  
                            idx = [id for id, x in enumerate(var_up_pos) if x in c and x is not var_up_pos[j]]                                              
                            in_cluster = [x for x in idx if distances[j][x] <= cutoff]
                            if set(in_cluster) == set(idx):
                               clusters_up[cl].append(var_up_pos[j])  
                               clusters_pdb[cl].append(var_pdb_pos[j])
                               clusters_wap[cl] += waps[i][j]
                               added = 1                                                                    
                        elif var_up_pos[j] in c:		
                            idx = [id for id, x in enumerate(var_up_pos) if x in c and x is not var_up_pos[i]]
                            in_cluster = [x for x in idx if distances[i][x] <= cutoff]  
                            if set(in_cluster) == set(idx):
                               clusters_up[cl].append(var_up_pos[i])    
                               clusters_pdb[cl].append(var_pdb_pos[i])
                               clusters_wap[cl] += waps[i][j]
                               added = 1
                    if not added:
                         clusters_up.append([var_up_pos[i], var_up_pos[j]])
                         clusters_pdb.append([var_pdb_pos[i], var_pdb_pos[j]])              
                         clusters_wap.append(waps[i][j])

   # Get out PyMol
   # pymol.cmd.quit()

   return (clusters_up, clusters_pdb, clusters_wap)
