library(TCGAbiolinks)
library(dplyr)
library(tidyverse)
library(readxl)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)

possibleEdges <- read.table("miRTarBaseEdges.txt",sep = "\n")$V1

raw.counts.mRNA <- read.csv(paste("normalized_mRNA.tsv",sep=""),sep="\t")

raw.counts.miRNA <- read.csv(paste("normalized_miRNA.tsv",sep=""),sep="\t")

cols <- colnames(raw.counts.mRNA)
cutoffPercentile <- 0.5
geneIds <- rownames(raw.counts.mRNA)
for(c in cols){
  vals <- raw.counts.mRNA[[c]]
  cutoff <- quantile(vals[!is.na(vals)],cutoffPercentile)
  sampleGenes <- filter(raw.counts.mRNA,raw.counts.mRNA[[c]]>cutoff)
  geneIds <- intersect(geneIds,rownames(sampleGenes))
}

cols <- colnames(raw.counts.miRNA)
miRNames <- rownames(raw.counts.miRNA)
for(c in cols){
  vals = raw.counts.miRNA[[c]]
  cutoff <- quantile(vals[!is.na(vals)],cutoffPercentile)
  sampleMiRs <- filter(raw.counts.miRNA,raw.counts.miRNA[[c]]>cutoff)
  miRNames <- intersect(miRNames,rownames(sampleMiRs))
}

print("Number of Genes")
print(length(geneIds))
print("Number of miRs")
print(length(miRNames))

samples <- intersect(colnames(raw.counts.mRNA),colnames(raw.counts.miRNA))
geneSampleDF <- dplyr::select(raw.counts.mRNA,all_of(samples))
geneSampleDF <- geneSampleDF %>% na.omit()
miRSampleDF <- dplyr::select(raw.counts.miRNA,all_of(samples))

keep <- c()
for(miR in miRNames){
  miRVals <- as.numeric(filter(miRSampleDF,rownames(miRSampleDF) == miR)[1,-1])
  if(length(miRVals) > 50){
    keep <- c(keep,miR)
  }
}
miRNames <- keep

corDF <- cor(t(geneSampleDF)[,geneIds],t(miRSampleDF)[,miRNames],method="spearman")
names <- expand.grid(miRNames,geneIds)
names <- unite(names,edges,sep=" ")
edges <- names[["edges"]]
correlations <- as.vector(t(corDF))

geneNames <- c()
for(gene in geneIds){
  if(as.character(gene) %in% keys(org.Hs.eg.db,keytype="ENTREZID")){
    geneName <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db,keys=c(as.character(gene)),column="SYMBOL",keytype="ENTREZID")[[1]])
  }
  else{
    geneName <- as.character(gene)
  }
  geneNames <- c(geneNames,geneName)
}

names <- expand.grid(miRNames,geneNames)
names <- unite(names,edges,sep=" ")
edges.wNames <- names[["edges"]]

corrCutoff <- 0.5
corrDF <- data.frame(edges,correlations)
corrDF.wNames <- data.frame(edges.wNames,correlations)
write.csv(corrDF,paste("BRCA_miRNA_mRNA_correlations.csv",sep=""))

corrDF <- corrDF %>% filter(correlations < 0)
vals <- corrDF[,2]
cutoff <- quantile(vals,corrCutoff)
edges <- filter(corrDF,correlations<cutoff)$'edges'
edges.wNames <- filter(corrDF.wNames,correlations<cutoff)$'edges.wNames'

realEdges2 <- c()
realEdges.wNames <- c()
for(i in 1:length(edges)){
  e <- edges[i]
  e.name <- edges.wNames[i]
  if(e %in% possibleEdges){
    realEdges2 <- c(realEdges2,e)
    realEdges.wNames <- c(realEdges.wNames,e.name)
  }
}

edgeList <- file(paste("BRCA_EdgeList.txt",sep=""))
writeLines(realEdges,edgeList)
close(edgeList)

edgeList <- file(paste("BRCA_EdgeList_wNames.txt",sep=""))
writeLines(realEdges.wNames,edgeList)
close(edgeList)
