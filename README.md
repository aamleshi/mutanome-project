# mutanome-project
Tool to identify mutation clusters in the 3D structures of proteins from MAF files

**D** : Directory
**[~]**: Data File


# Pre-Process:
- \> `DownloadMaf.py` (Download Maff files / vertify their existance) 
  - **D** Maf_Files
    - **D** Cancer type 1  
        - **[~]** pipeline 1
        - **[~]** pipeline 2 
        - **[~]** pipeline 3 
        - **[~]** pipeline 4 
    - **D** Cancer type n
        - **[~]** pipeline 1
- \> `GroupGenes.py` (Create a set of dataframes for each gene for each pipeline) (TODO)
  - **[~]** gene_table.pkl
  - **D** Pipeline1
    - **D** Chr 1
      - **[~]** geneA.pkl
      - **[~]** geneB.pkl
- \> `Download_PDB` (TODO)
  - **D** PDB_Files
    - **D** Chr 1
      - **D** GeneA
        - **[~]** Residue1.pdb
        - **[~]** Residue2.pdb
  # Analysis:
  TODO
  
