# Construction and Analysis of bipartite miRNA-mRNA network
## Prerequisites
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
- R 4.3.3
  - tidyverse 2.0.0
  - AnnotationDbi 1.64.1
  - Biocondutor 3.21
  - org.Hs.eg.db 3.18

To run the code, the following files are also needed:
- List of experimentally validated interactions from miRTarBase
  - Provided as mirRTarBaseEdges.txt
- Normalized miRNA and mRNA expression from TCGA
  - Available through NIH Genomic Data Commons TCGA Pan-Cancer Atlas (Hoadley, K. A. et al. (2018) https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin)
- GMT files of pathways of interest
  - KEGG hsa Pathways
  - Biocarta Pathways
  - DAVID Bioinformatics Functional Annotation Pathways
  - mSigDB Pathways 

## Network Construction Pipeline

## Network Community Analysis Pipeline
