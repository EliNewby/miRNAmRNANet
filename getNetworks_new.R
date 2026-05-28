library(TCGAbiolinks)
library(dplyr)
library(tidyverse)
library(readxl)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ppcor)

#miRTarBaseData = read_excel("Documents/mRNAmiRNANets/hsa_MTI.xlsx")
#ids <- as.vector(miRTarBaseData %>% dplyr::select(`Target Gene (Entrez ID)`))
#possibleEdges <- c()
#for(r in row.names(miRTarBaseData)[1:100]){
#  miR <- miRTarBaseData[r,"miRNA"]
#  gene <- miRTarBaseData[r,"Target Gene (Entrez ID)"]
#  if(!(paste(miR,gene) %in% possibleEdges)){
#    possibleEdges <- c(possibleEdges,paste(miR,gene))
#  }
#}
possibleEdges <- read.table("Documents/mRNAmiRNANets/possibleEdges.txt",sep = "\n")$V1

#library(TCGAbiolinks)
#library(dplyr)
#library(tidyverse)
#library(readxl)
#library(biomaRt)
#library(AnnotationDbi)
#library(org.Hs.eg.db)
#library(ppcor)
#query.mRNA <- GDCquery(project = "TCGA-LAML", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
#GDCdownload(query.mRNA)
#raw.counts.mRNA <- GDCprepare(query = query.mRNA, summarizedExperiment = FALSE)

#raw.counts.miRNA <- read_csv("Documents/mRNAmiRNANets/BRCA_cleaned_miRNA_mature.csv")
#query.miRNA <- GDCquery(project = "TCGA-LAML", data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification")
#GDCdownload(query.miRNA)
#raw.counts.miRNA <- GDCprepare(query = query.miRNA, summarizedExperiment = FALSE)

#cancers <- c("LUAD","COAD","OV","BLCA","LGG","KIRC","LUSC","ACC","CESC","LAML")#c("BRCA","LUAD","COAD","OV","BLCA","LGG","KIRC","LUSC","LAML","ACC","CESC","CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRP","LIHC","MESO","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
#cancers <- c("CHOL","DLBC","ESCA","GBM","HNSC","KICH","KIRP","LIHC","MESO","PAAD","PCPG")
#cancers <- c("PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
cancers <- c("BRCA")

for(cancer in cancers){
  print(cancer)
  raw.counts.mRNA <- read.csv(paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/mRNA/normalized_mRNA.tsv",sep=""),sep="\t")
  cols <- colnames(raw.counts.mRNA)
  if(length(cols) < 50){
    next
  }
  #lastColName <- cols[length(cols)]
  #lastCol <- raw.counts.mRNA[[lastColName]]
  #v1 <- c()
  #v2 <- c()
  #for(x in lastCol){
  #  vals <- str_split(x,"\t")[[1]]
  #  v1 <- c(v1,as.numeric(vals[1]))
  #  v2 <- c(v2,as.numeric(vals[2]))
  #}
  #raw.counts.mRNA[[lastColName]] <- NULL
  #raw.counts.mRNA$v1 <- v1
  #raw.counts.mRNA$v2 <- v2
  #colnames(raw.counts.mRNA) <- c("gene_id",cols)
  
  raw.counts.miRNA <- read.csv(paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/miRNA/normalized_miRNA.tsv",sep=""),sep="\t")
  cols <- colnames(raw.counts.miRNA)
  #lastColName <- cols[length(cols)]
  #lastCol <- raw.counts.miRNA[[lastColName]]
  #v1 <- c()
  #v2 <- c()
  #for(x in lastCol){
  #  vals <- str_split(x,"\t")[[1]]
  #  v1 <- c(v1,as.numeric(vals[1]))
  #  v2 <- c(v2,as.numeric(vals[2]))
  #}
  #raw.counts.miRNA[[lastColName]] <- NULL
  #raw.counts.miRNA$v1 <- v1
  #raw.counts.miRNA$v2 <- v2
  #colnames(raw.counts.miRNA) <- c("mirna_id",cols)
  
  cols <- colnames(raw.counts.mRNA)
  cutoffPercentile <- 0.5
  geneIds <- rownames(raw.counts.mRNA)
  #oldGeneIds <- rownames(raw.counts.mRNA)
  quantDF <- data.frame(geneIds)
  i = 1
  for(c in cols){
    if(i%%100 == 0){
      print(i/length(cols)*100)
      #print(length(geneIds))
    }
    i <- i+1
    vals <- raw.counts.mRNA[[c]]
    quants <- ecdf(vals)(vals)
    quantDF[[c]] <- quants
    cutoff <- quantile(vals[!is.na(vals)],cutoffPercentile)
    sampleGenes <- filter(raw.counts.mRNA,raw.counts.mRNA[[c]]>cutoff)
    geneIds <- intersect(geneIds,rownames(sampleGenes))
  }
  
  minQuants <-c()
  nearMinQuants <- c()
  for(gene in rownames(raw.counts.mRNA)){
    vals <- as.numeric(filter(quantDF,geneIds==gene)[1,-1])
    min <-  min(vals)
    minQuants <- c(minQuants,min)
    nearMin <- quantile(vals,0.1,na.rm=T)[[1]]
    nearMinQuants <- c(nearMinQuants,nearMin)
  }
  minQuantDF <- data.frame(rownames(raw.counts.mRNA),minQuants,nearMinQuants)
  #geneIds <- (minQuantDF %>% filter(minQuants >= cutoffPercentile))[,1]
  #geneIds <- (minQuantDF %>% filter(nearMinQuants > cutoffPercentile))[,1]
  #write.csv(quantDF,paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/mRNA/AllQuantiles.csv",sep=""))
  
  cols <- colnames(raw.counts.miRNA)
  miRNames <- rownames(raw.counts.miRNA)
  quantDF <- data.frame(miRNames)
  
  i = 1
  for(c in cols){
    if(i%%100 == 0){
      print(i/length(cols)*100)
      #print(length(miRNames))
      }
    i <- i+1
    vals = raw.counts.miRNA[[c]]
    quants <- ecdf(vals)(vals)
    quantDF[[c]] <- quants
    cutoff <- quantile(vals[!is.na(vals)],cutoffPercentile)
    sampleMiRs <- filter(raw.counts.miRNA,raw.counts.miRNA[[c]]>cutoff)
    miRNames <- intersect(miRNames,rownames(sampleMiRs))
  }
  
  minQuants <- c()
  nearMinQuants <- c()
  miRs <- c()
  for(miR in rownames(raw.counts.miRNA)){
    vals <- as.numeric(filter(quantDF,miRNames==miR)[1,-1])
    min <- min(vals,na.rm=TRUE)
    if(!is.infinite(min)){
      minQuants <- c(minQuants,min)
      nearMin <-quantile(vals,0.1)[[1]]
      nearMinQuants <- c(nearMinQuants,nearMin)
      miRs <- c(miRs,miR)
    }
  }
  minQuantDF <- data.frame(miRs,minQuants,nearMinQuants)
  #miRNames <- (minQuantDF %>% filter(minQuants >= cutoffPercentile))[,1]
  #miRNames <- (minQuantDF %>% filter(nearMinQuants > cutoffPercentile))[,1]
  #write.csv(quantDF,paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/miRNA/AllQuantiles.csv",sep=""))
  
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
  
  #geneNames <- c()
  #for(gene in geneIds){
  #  if(as.character(gene) %in% keys(org.Hs.eg.db,keytype="ENTREZID")){
  #    geneName <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db,keys=c(as.character(gene)),column="SYMBOL",keytype="ENTREZID")[[1]])
  #  }
  #  else{
  #    geneName <- as.character(gene)
  #  }
  #  geneNames <- c(geneNames,geneName)
  #}

  #names <- expand.grid(miRNames,geneNames)
  #names <- unite(names,edges,sep=" ")
  #edges.wNames <- names[["edges"]]
  edges <- c()
  #edges.wNames <- c()
  correlations <- c()
  #pVals <- c()
  i <- 1
  for(gene in geneIds){
    if(i%%10 == 0){
      print(i/length(geneIds)*100)
    }
    i <- i+1
    geneVals <- as.numeric(filter(geneSampleDF,rownames(geneSampleDF) == gene))
    for(miR in miRNames){
      miRVals <- as.numeric(filter(miRSampleDF,rownames(miRSampleDF) == miR))
      geneVals <- geneVals[!is.na(miRVals)]
      miRVals <- miRVals[!is.na(miRVals)]
      if(length(miRVals) > 50){
        correlation <- cor.test(miRVals,geneVals,method="spearman")
        pVal <- correlation$p.value
        rho <- as.numeric(correlation$estimate)
        correlations <- c(correlations,rho)
        pVals <- c(pVals,pVal)
        edges <- c(edges,paste(miR,gene))
        #if(as.character(gene) %in% keys(org.Hs.eg.db,keytype="ENTREZID")){
        #  geneName <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db,keys=c(as.character(gene)),column="SYMBOL",keytype="ENTREZID")[[1]])
        #}
        #else{
        #  geneName <- as.character(gene)
        #}
        #edges.wNames <- c(edges.wNames,paste(miR,geneName))
      }
    }  
  }
  
  corrCutoff <- 0.5
  #corrDF <- data.frame(edges,correlations,pVals)
  corrDF <- data.frame(edges,correlations)
  #corrDF.wNames <- data.frame(edges.wNames,correlations,pVals)
  #write.csv(corrDF,paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/miRNA_mRNA_correlations_nearMin.csv",sep=""))
  
  corrDF <- corrDF %>% filter(correlations < 0)
  vals <- corrDF[,2]
  cutoff <- quantile(vals,corrCutoff)
  edges <- filter(corrDF,correlations<cutoff)$'edges'
  #edges.wNames <- filter(corrDF.wNames,correlations<cutoff)$'edges.wNames'
  
  realEdges2 <- c()
  #realEdges.wNames <- c()
  for(i in 1:length(edges)){
    if(i%%1000 == 0){
      print(i/length(edges)*100)
    }
    e <- edges[i]
    #e.name <- edges.wNames[i]
    if(e %in% possibleEdges){
      realEdges2 <- c(realEdges2,e)
      #realEdges.wNames <- c(realEdges.wNames,e.name)
    }
  }
  
  #edgeList <- file(paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/",cancer,"_EdgeList_matrix.txt",sep=""))
  #writeLines(realEdges,edgeList)
  #close(edgeList)
  
  #edgeList <- file(paste("Documents/mRNAmiRNANets/tcga-data-2-2025/",cancer,"/",cancer,"_EdgeList_nearMin_wNames.txt",sep=""))
  #writeLines(realEdges.wNames,edgeList)
  #close(edgeList)
}