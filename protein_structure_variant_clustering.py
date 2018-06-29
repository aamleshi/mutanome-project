
#!/usr/bin/python

import subprocess
import os
import sys
import glob
import re
import urllib
import pickle
from Bio import SwissProt
from classes.protein_info import PDB

from classes.protein_info import ProtVar
from variants_functions import check_pdb_coverage
from variants_functions import create_clusters
from variants_functions import read_pdb

# Print USAGE
if len(sys.argv) < 4:
   print """\
   There are some missing arguments.
   Usage: python maf_variants_protein_structures_clustering_integrated.py PATH OUT_PATH N_SAMPLES PROTEINS
        PATH: Path where MAF files are included   
        OUT_PATH: Path where partial and final results will be stored
        N_SAMPLES: Number of samples considered for recurrent variants
        PROTEINS: File with list of proteins we are interested in
   """
   sys.exit()  
else:
   
   leng = len(sys.argv)
   mafpaths = glob.glob('*' + sys.argv[1])
   print mafpaths
   outpath = sys.argv[2]
   nvar = int(sys.argv[3])
   protfile = sys.argv[4]
   ncallers = 1
   pdbsource = "ASSEMBLY"

# Read file or list of proteins
if '.' in protfile:
   with open(protfile, 'r') as f:
      protids = f.read().splitlines()
elif len(protfile) > 0: 
   protids = protfile.split(',') 
else:
   protids = ""

###################################################################################
#
# STEP 1: Reading MAFs, integrating and filtering variants
#
###################################################################################
obj_file = os.path.basename(protfile).replace('.txt', '.pckl')
print obj_file

# Read file object if exists
if os.path.exists(obj_file):
   with open(obj_file) as fobj:
      variants_list = pickle.load(fobj)
# Calculate if not
else:
   variants_list = {}
   for caller in mafpaths:
      print ('Finding recurrent variants from MAFs files in ' + caller + '...')
      mafs = sorted(glob.glob(caller + '/*.maf'))
      for file in mafs:

         # Read cancer type and caller
         cancer = file.split('.')[1]
         pipeline = file.split('.')[2]

         # Create a list of cases
         cases_list = []

         # Reading each maf file
         with open(file) as f:
            for line in f:

               # Read columns for each variant in the MAF file
               columns = line.split('\t')

               # Filter empty rows and headers
               if len(columns)>2 and columns[0] != "Hugo_Symbol":               

                  # Save case if not considered before            
                  case = columns[15] + ' ' + columns[16]
                  if case not in cases_list:
                      cases_list += [case]

                  # Consider variants only in the protein list
                  if not protids or columns[67] in protids:
    
                     # STEP 1: Filtering variants
                     # 1) Missense mutations (SNPs)
                     # 2) Filter "PASS"
                     # 3) At least nvar recurrences for each mutation at protein level
                     # 4) Only "Somatic" mutations (just in case)             
                     # 5) Include only non empty Uniprot  ID
                     if columns[8] == "Missense_Mutation" and columns[9] == "SNP" and columns[25] == "Somatic" and (columns[108] == "PASS" or columns[108] == "common_variant") and columns[67] != "" :
                   
                         # Retrieve location of variant in protein
                         location = columns[67] + ":" + columns[53] 

                         # Format variant position and AAs
                         posvar = columns[53].split('/')
                         formatted_var = columns[54].replace('/',posvar[0])                        

                         if location in variants_list:
                         
                             # Check cases
                             if case not in variants_list[location][4]:
                                 #if formatted_var not in variants_list[location][3]:
                                 #     variants_list[location][3] += [formatted_var]
                                 aa_vars = variants_list[location][3].split("/")
                                 aa_vars[0] = aa_vars[0][-1]
                                 if formatted_var[-1] not in aa_vars:
                                      variants_list[location][3] += "/" + formatted_var[-1] 
                                 variants_list[location][4] += [case]
                                 variants_list[location][5] += 1

                                 # Check projects
                                 if cancer in variants_list[location][8]:
                                      cancer_idx = variants_list[location][8].index(cancer)
                                      variants_list[location][9][cancer_idx] += 1
                                 else:
                                      variants_list[location][8] += [cancer]
                                      variants_list[location][9] += [1]


                             # Check pipelines
                             if pipeline not in variants_list[location][6]:
                                 variants_list[location][6] += [pipeline]
                                 variants_list[location][7] += 1 

                         else:                   
                             variants_list[location] = [columns[67], columns[60], columns[53], formatted_var, [case], 1, [pipeline], 1, [cancer], [1]]           

         # Close caller file
         f.close()  

         # Print number of cases
         print "CASES:\t{0}\t{1}\t{2}".format(caller, cancer, len(cases_list)) 

   # Filtering variants which are in less than 2 callers 
   variants_list = { key:value for key, value in variants_list.items() if value[5] >= nvar and value[7] >= ncallers}

   # Save variant_list for future use
   fvar = open(obj_file, 'wb')
   pickle.dump(variants_list, fvar)
   fvar.close()

#######################################################################################
#
# STEP 2: Retrieving SWISSPROT associated to missense variants in the filtered variants
#
#######################################################################################
print ('Retrieving protein information...')
protvar_list = []
#fstats = open(outpath + 'variants_stats.txt', 'w')
#fstats.write("SELECTION\tUP_ID\tGENE\tVAR_POS\tPROJECTS\tFREQUENCY\n")
for key, variant in variants_list.iteritems():   
   
   # Consider only if Uniprot ID is not empty
   if variant[0] != "":
       
       # Reading protein in Uniprot
       if ',' in variant[0]:
           variant[0] = variant[0][:variant[0].index(',')]       
 
       # Downloading Uniprot file if doesn't exist
       prot_file = "uniprot/" + variant[0] + ".txt"
       if not os.path.exists(prot_file):
           print ("   Downloading " + prot_file + " ...")
           data = urllib.urlopen("http://www.uniprot.org/" + prot_file)
           fw = open(prot_file, 'w+')
           fw.write(data.read()) 
           fw.close()
         
       # If file is empty, save it and continue to next iteration
       if os.stat(prot_file).st_size == 0:      
           # fstats.write("EMPTY_FILE\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0],variant[1],key,variant[8],variant[5]))     
           continue

       # Retrieve information from Uniprot: sequence length and AA in variant
       handle = open(prot_file)
       uniprot = SwissProt.parse(handle).next()
       length = uniprot.sequence_length

       # Reading positions and total length
       posvar=variant[2].split("/")

       # Check if AA annotation and sequence length match (same isoform)                
       if length == int(posvar[1]):
            aa = uniprot.sequence[int(posvar[0])-1]
            if aa == variant[3][0][0]:   # CHECK THIS TO CONSIDER ALL POSSIBILITIES              	
                idx =  [i for i,x in enumerate(protvar_list) if x.protid == variant[0]]               
                chrpos = key.split(':') 
                formproj = ["{}:{}".format(a_, b_) for a_, b_ in zip(variant[8], variant[9])]    
                if not idx:
                    proteinInfo = ProtVar(variant[0], chrpos[0], [chrpos[1]], "missense_mutation", [variant[5]], variant[1], [int(posvar[0])], [int(posvar[1])], [variant[3]], [formproj])                
                    protvar_list.append(proteinInfo)
                    # fstats.write("SELECT\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0],variant[1], key, variant[8], variant[5]))     
                elif chrpos[1] not in protvar_list[idx[0]].pos:		                      
                    protvar_list[idx[0]].pos += [chrpos[1]]
                    protvar_list[idx[0]].freq += [variant[5]]
                    protvar_list[idx[0]].aa_pos += [int(posvar[0])]                      
                    protvar_list[idx[0]].aa_var += [variant[3]]     
                    protvar_list[idx[0]].project += [formproj]               		                
                    # fstats.write("SELECT\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0],variant[1], key, variant[8], variant[5]))
                    # print record.CHROM + ":{0}".format(record.POS) + "\t" +  info_split[1] + "\t" + info_split[3] + "\t" + info_split[14] + "\t" + info_split[15] + "\t" + info_split[30] + "\t" + "\t" + info_split[33]
                # else:
                    # fstats.write("ALREADY_ANN\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0],variant[1],key,variant[8],variant[5]))
            # else:
                # fstats.write("DIF_AA_ISOFORM\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0], variant[1], key,variant[8],variant[5]))
       # else:
           # fstats.write("DIF_LEN_ISOFORM\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0],variant[1],key,variant[8],variant[5]))

       # Close Uniprot file
       handle.close()

   # else:
       # fstats.write("NO_UPID\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(variant[0],variant[1],key, variant[8],variant[5]))

# Close variant statistics file
# fstats.close()


###################################################################################
#
# STEP 3: Determining PDB structures covering proteins and their variants
#
###################################################################################
print "Determining PDB structures covering variants..."
fprot = open(outpath + 'proteins_stats.txt', 'w')
fprot.write("PROTEIN\tGENE\tNVAR\tVAR_POS\tVAR_AA\tFREQUENCY\tPROJECT\n")
fncov = open(outpath + 'non_covered_variants_stats.txt', 'w')
fncov.write("PROTEIN\tGENE\tREASON\tNVAR\tEXC_NVAR\tEXC_VAR_POS\n")
fpdb = open(outpath + 'pdbs_stats.txt', 'w')
fpdb.write("PDB\tNVAR\tVAR_POS_UP\tVAR_POS_PDB\tVAR_AA\tPROJECT\tFREQUENCY\tPROTEIN\tGENE\n")
selected_pdbs = []
for protein in protvar_list:

   # For statistics purpose
   fprot.write("{0}\t{1}\t{2}\t\"{3}\"\t\"{4}\"\t\"{5}\"\t\"{6}\"\n".format(protein.protid, protein.gene, len(protein.aa_pos), protein.aa_pos, protein.aa_var, protein.freq, protein.project))

   # Reading protein in Uniprot
   prot_file = "uniprot/" + protein.protid + ".txt"
   if os.path.exists(prot_file):
      handle = open(prot_file)

   # Finding PDB records in Uniprot file
   record = SwissProt.parse(handle).next()
   UPseq = record.sequence
   variantsPos = protein.aa_pos
   variants = protein.aa_var
   projects = protein.project
   frequencies = protein.freq
   pdb_list = []
   for feature in record.cross_references:
       if feature[0] == "PDB":
           seqpos = feature[4].split(",")
           if feature[4] == "-":
               continue
           m = re.match(r"(.+)=(\d+)-(\d+)",seqpos[0])
           if hasattr(m, 'group'):
              chain = m.group(1)
              init = int(m.group(2))
              ending = int(m.group(3))
              count = 0
              incvar = []
              incvarpos = []
              incproj = []
              incfreq = []
              for pos in variantsPos:
                  if pos >= init and pos <= ending and pos not in incvarpos:
                      count+=1
                      incvarpos = incvarpos + [pos]
                      varposit = variantsPos.index(pos)
                      incvar = incvar + [variants[varposit]]
                      incproj = incproj + [projects[varposit]]
                      incfreq = incfreq + [frequencies[varposit]]
              length = ending - init + 1
              pdb_list.append(PDB(feature[1], chain, init, ending, length, count, incvar, incvarpos, incvarpos, protein.protid, protein.gene, incfreq, incproj))  

   # Select those PDB maximaxing the sequence and variants coverage
   if len(pdb_list) > 0:
      sorted_pdbs = sorted(pdb_list, key=lambda pdb: (pdb.num_variants, pdb.length), reverse=True)
      remainVar = set(variantsPos)
      nvar_pdb_min = 2

      # Select PDB if:
      # 1. Include variants not included before
      # 2. Have a identificable sequence (FASTA file)
      # 3. Include at least 2 variants in total
      for pdb in sorted_pdbs:
           included = remainVar & set(pdb.var_up_pos)
           if included and len(variantsPos) >= nvar_pdb_min:

              # Check effective coverage of variants in this PDB structure       
              structures = read_pdb(pdb.name, pdbsource)
              max_num_var = 0
              updated_pdb = []
              if structures:
                 for struct in structures:
                    mod_pdb = check_pdb_coverage(pdb, UPseq, struct, outpath, pdbsource)
                    if mod_pdb:
                       st_num_var = len(mod_pdb.var_pdb_pos)
                       #if st_num_var > max_num_var:
                       #     max_num_var = st_num_var
                       #     updated_pdb = mod_pdb                            
                       selected_pdbs.append(mod_pdb)
                       fpdb.write("{0}\t{1}\t\"{2}\"\t\"{3}\"\t\"{4}\"\t\"{5}\"\t\"{6}\"\t{7}\t{8}\n".format(mod_pdb.name, len(mod_pdb.var_pdb_pos), mod_pdb.var_up_pos, mod_pdb.var_pdb_pos, mod_pdb.variants, mod_pdb.project, mod_pdb.freq, mod_pdb.protein, mod_pdb.gene))
        
              #if updated_pdb:
              #    selected_pdbs.append(updated_pdb)
              #    fpdb.write("{0}\t{1}\t\"{2}\"\t\"{3}\"\t\"{4}\"\t\"{5}\"\t\"{6}\"\t{7}\t{8}\n".format(updated_pdb.name, len(updated_pdb.var_pdb_pos), updated_pdb.var_up_pos, updated_pdb.var_pdb_pos, updated_pdb.variants, updated_pdb.project, updated_pdb.freq, updated_pdb.protein, updated_pdb.gene))
              #    remainVar = remainVar - set(updated_pdb.var_up_pos)

              #if not remainVar:
              #    break

      # Variants not covered by any PDB
      if len(remainVar) > 0:
          fncov.write('{0}\t{1}\tPDB_NOT_COVERING\t{2}\t{3}\t\"{4}\"\n'.format(protein.protid, protein.gene, len(set(variantsPos)), len(remainVar), list(remainVar)))

   # Proteins without any PDB associated
   else:
      fncov.write('{0}\t{1}\tNOT_PDB\t{2}\t{3}\t\"{4}\"\n'.format(protein.protid, protein.gene, len(protein.aa_pos), len(protein.aa_pos), protein.aa_pos))

   # Close Uniprot fie
   handle.close()

# Close statistics files
fprot.close()
fncov.close()
fpdb.close()


###################################################################################
#
# STEP 4: Calculating distances of variants in PDB and find clustering
#
###################################################################################   
del protvar_list
sorted_selected_pdbs = sorted(selected_pdbs, key=lambda pdb: pdb.num_variants, reverse=True)

# Calculate distance arrays
print ("Calculating distances and variant clustering in 3D structures...")
fdis = open(outpath + 'variant_clustering.txt', 'w')
fdis.write("PDB\tPROTEIN\tGENE\tNVAR_CLUSTER\tCLUSTER_UP_POS\tCLUSTER_PDB_POS\tFREQUENCIES\tPROJECTS\tWAP_SCORE\tCHAIN\n")
cutoff = 14
for pdb in sorted_selected_pdbs:
   if (pdb.num_variants>1):

      print "      Looking for clusters in " + pdb.name + " ..."
      clustersUP, clustersPDB, clustersWAP = create_clusters(pdb, cutoff)

      for c in range(0, len(clustersUP)):
         clusterup = clustersUP[c]
         formatted_cluster = [pdb.variants[idx] for idx, x in enumerate(pdb.var_up_pos) if x in clusterup]
         projects_cluster = [pdb.project[idx] for idx, x in enumerate(pdb.var_up_pos) if x in clusterup]
         freq_cluster = [pdb.freq[idx] for idx, x in enumerate(pdb.var_up_pos) if x in clusterup]
         fdis.write("{0}\t{1}\t{2}\t{3}\t\"{4}\"\t\"{5}\"\t\"{6}\"\t\"{7}\"\t{8}\t{9}\n".format(pdb.name, pdb.protein, pdb.gene, len(clustersUP[c]), formatted_cluster, clustersPDB[c], freq_cluster, projects_cluster, clustersWAP[c], pdb.chain))

fdis.close()
