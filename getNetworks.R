library(TCGAbiolinks)
library(dplyr)
library(tidyverse)
library(readxl)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rstatix)
library(infotheo)
library(igraph)

#miRTarBaseData = read_excel("hsa_MTI.xlsx")
#ids <- as.vector(miRTarBaseData %>% dplyr::select(`Target Gene (Entrez ID)`))
#possibleEdges <- c()
#for(r in row.names(miRTarBaseData)){
#  miR <- miRTarBaseData[r,"miRNA"]
#  gene <- miRTarBaseData[r,"Target Gene (Entrez ID)"]
#  if(!(paste(miR,gene) %in% possibleEdges)){
#    possibleEdges <- c(possibleEdges,paste(miR,gene))
#  }
#}
possibleEdges <- read.table(" possibleEdges.txt",sep = "\n")$V1
cancer <- "BRCA"
subtypeData <- TCGAquery_subtype("BRCA")
subtypeData <- subtypeData %>% mutate(patient = str_replace_all(patient, "-", "."))

raw.counts.mRNA <- read.csv(paste("BRCA_normalized_mRNA.tsv",sep=""),sep="\t")
cols <- colnames(raw.counts.mRNA)
colKeep <- c()
for(c in cols){
  temp <- str_split(c,"\\.")
  if(startsWith(temp[[1]][4],"01")){
    colKeep <- c(colKeep,TRUE)
  }
  else{colKeep <- c(colKeep,FALSE)}
}
cols <- cols[colKeep]
raw.counts.mRNA <- raw.counts.mRNA %>% dplyr::select(all_of(cols))
colNames <- lapply(cols,str_sub,start=1,end=12)
colnames(raw.counts.mRNA) <- colNames

raw.counts.miRNA <- read.csv(paste("BRCA_normalized_miRNA.tsv",sep=""),sep="\t")
cols <- colnames(raw.counts.miRNA)
colKeep <- c()
for(c in cols){
  temp <- str_split(c,"\\.")
  if(startsWith(temp[[1]][4],"01")){
    colKeep <- c(colKeep,TRUE)
  }
  else{colKeep <- c(colKeep,FALSE)}
}
cols <- cols[colKeep]
raw.counts.miRNA <- raw.counts.miRNA %>% dplyr::select(all_of(cols))
colNames <- lapply(cols,str_sub,start=1,end=12)
colnames(raw.counts.miRNA) <-  colNames

samples <- intersect(colnames(raw.counts.mRNA),colnames(raw.counts.miRNA))
raw.counts.mRNA <- dplyr::select(raw.counts.mRNA,all_of(samples))
raw.counts.mRNA <- raw.counts.mRNA %>% na.omit()
raw.counts.miRNA <- dplyr::select(raw.counts.miRNA,all_of(samples))
raw.counts.miRNA <- raw.counts.miRNA %>% na.omit()

### Get Genes/MiRNAs with expression above Q3

quantileGeneExp <- apply(raw.counts.mRNA,1,quantile,0.75)
geneIds <- names(quantileGeneExp[quantileGeneExp>quantile(quantileGeneExp,0.75)])

write_lines(geneIds,paste("BRCA_Genes.txt",sep=""))

quantileMiRExp <- apply(raw.counts.miRNA,1,quantile,0.75)
miRNames <- names(quantileMiRExp[quantileMiRExp>quantile(quantileMiRExp,0.75)])

names <- expand.grid(miRNames,geneIds)
names <- unite(names,edges,sep=" ")
edgeNames <- names[["edges"]]

geneDF <- data.frame(t(raw.counts.mRNA)[,geneIds])
miRDF <- data.frame(t(raw.counts.miRNA)[,miRNames])
mergedDF <- data.frame(geneDF,miRDF)
geneDFDiscrete <- discretize(geneDF)
miRDFDiscrete <- discretize(miRDF)

vars1 <- colnames(geneDF)
vars2 <- colnames(miRDF)

corDF <- cor(geneDF,miRDF,method="spearman")
correlations <- as.vector(t(corDF))
corrDF <- data.frame(edgeNames,correlations)
write.csv(corrDF,paste("BRCA_miRNA_mRNA_correlations.csv",sep=""))

thresholds <- seq(-0.3,-0.01,by=0.01)
for (t in thresholds){
  print(t)
  temp <- filter(corrDF,correlations<t)
  edges <- temp$edgeNames
  realEdges <- intersect(edges,possibleEdges)
  if(length(realEdges)>0){
    edgeDF <- read.table(text=realEdges,col.names=c("miR","mRNA"))
    G <- graph_from_data_frame(edgeDF,directed=F)
    GLCC <- largest_component(G)
    if(vcount(G)>50){
      fracOfNodes <- vcount(GLCC)/vcount(G)
    }
    else{
      fracOfNodes <- 0
    }
    if(fracOfNodes == 1){
      thresh <- t
      break
    }
  }
}
edgeList <- file(paste("BRCA_EdgeList.txt",sep=""))
writeLines(realEdges,edgeList)
close(edgeList)


### Repeat for BRCA Subtypes

for(subtype in c("Basal","Her2","LumA","LumB")){
  print(subtype)
  subtypeSamples <- (subtypeData %>% filter(BRCA_Subtype_PAM50 == subtype))$patient
  
  geneExpressionDF <- raw.counts.mRNA %>% select(intersect(subtypeSamples,samples))
  miRExpressionDF <- raw.counts.miRNA %>% select(intersect(subtypeSamples,samples))
  
  ### Get Genes/MiRNAs with expression above Q3
  
  quantileGeneExp <- apply(geneExpressionDF,1,quantile,0.75)
  geneIds <- names(quantileGeneExp[quantileGeneExp>quantile(quantileGeneExp,0.75)])
  
  write_lines(geneIds,paste("BRCA_",subtype,"_Genes.txt",sep=""))
  
  quantileMiRExp <- apply(miRExpressionDF,1,quantile,0.75)
  miRNames <- names(quantileMiRExp[quantileMiRExp>quantile(quantileMiRExp,0.75)])
  
  names <- expand.grid(miRNames,geneIds)
  names <- unite(names,edges,sep=" ")
  edgeNames <- names[["edges"]]
  
  geneDF <- data.frame(t(geneExpressionDF)[,geneIds])
  miRDF <- data.frame(t(miRExpressionDF)[,miRNames])
  mergedDF <- data.frame(geneDF,miRDF)
  geneDFDiscrete <- discretize(geneDF)
  miRDFDiscrete <- discretize(miRDF)
  
  vars1 <- colnames(geneDF)
  vars2 <- colnames(miRDF)
  
  corDF <- cor(geneDF,miRDF,method="spearman")
  correlations <- as.vector(t(corDF))
  corrDF <- data.frame(edgeNames,correlations)
  write.csv(corrDF,paste("BRCA_",subtype,"_miRNA_mRNA_correlations.csv",sep=""))
  
  thresholds <- seq(-0.3,-0.01,by=0.01)
  for (t in thresholds){
    print(t)
    temp <- filter(corrDF,correlations<t)
    edges <- temp$edgeNames
    realEdges <- intersect(edges,possibleEdges)
    if(length(realEdges)>0){
      edgeDF <- read.table(text=realEdges,col.names=c("miR","mRNA"))
      G <- graph_from_data_frame(edgeDF,directed=F)
      GLCC <- largest_component(G)
      if(vcount(G)>50){
        fracOfNodes <- vcount(GLCC)/vcount(G)
      }
      else{
        fracOfNodes <- 0
      }
      if(fracOfNodes == 1){
        thresh <- t
        break
      }
    }
  }
  edgeList <- file(paste("BRCA_",subtype,"_EdgeList.txt",sep=""))
  writeLines(realEdges,edgeList)
  close(edgeList)
  
}