# Construction and Analysis of bipartite miRNA-mRNA network
Python and R code for the creation and community-based analysis of bipartite miRNA-mRNA networks. 

Please cite
Multiomic network analysis reveals conserved and subtype-specific cooperative microRNA regulators in breast cancer
Eli Newby, Elijah Davis, Andrew Dhawan
bioRxiv 2025.08.07.669176; doi: https://doi.org/10.1101/2025.08.07.669176
  
## Network Construction Pipeline
The getNetworks.R file contstructs a bipartite miRNA-mRNA network based on mRNA and miRNA data from breast cancer patient data. It takes as input the miRNA and mRNA expression data (obtained from Hoadley, K. A. et al. Cell-of-Origin patterns dominate the molecular classification
of 10,000 tumors from 33 types of cancer. Cell 173, 291–304.e6 (2018)) as well as a list of experimentally verified interactions from miRTarBase (listed in possibleEdges.txt), and constructs a bipartite network based on highly negativley correlated miRNA and mRNAs. It outputs the network as an edgelist. The output edgelist files for breast cancer (BRCA) as well as 29 other cancer types are available through the NetworkEdgelists folder.

## Network Community Analysis Pipeline
The graphToFunctions.py file runs an analysis pipeline to determine functional annotations of a network. It takes as input an edgelist and analyzes the communities within this network for significant functional annotations. It outputs a file that lists the identified communities within the bipartite network, the functional annotations identified for each community, the proteomic associations of each community, and a summary file that clusters signficant annotations to determine overall functions withing each community of the network. 

The validateCommunities.py file runs a qstest (https://github.com/skojaku/qstest) on the generated network and identified community to determine which communities are statitically significant.

The synergyPPI.py file runs a combinatorics analysis to identify which miRNA combination is the most effective at driving an input biological pathway.

## Requirements
This code was written and compiled on:<br>
- Python 3.12.3
  -  https://github.com/genisott/pycondor
  -  https://github.com/skojaku/qstest
  -  rpy2 3.5.12
  -  networkx 3.3
  -  pandas 2.2.2
  -  numpy 1.26.4
  -  matplotlib 3.9.2
  -  scipy 1.13.1
  -  pycondor (https://github.com/genisott/pycondor)
- R 4.3.3
  - tidyverse 2.0.0
  - AnnotationDbi 1.64.1
  - Biocondutor 3.21
  - org.Hs.eg.db 3.18

To run the code, the following files are also needed:
- List of experimentally validated interactions from miRTarBase
  - Provided as possibleEdges.txt
- Normalized miRNA and mRNA expression from TCGA
  - Available through NIH Genomic Data Commons TCGA Pan-Cancer Atlas (Hoadley, K. A. et al. (2018) https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin)
- GMT files of pathways of interest from mSigDB Pathways
