if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("tidyverse",quietly=TRUE))
  BiocManager::install("tidyverse")
if(!require("fgsea",quietly=TRUE))
  BiocManager::install("fgsea")
if(!require("AnnotationDbi",quietly=TRUE))
  BiocManager::install("AnnotationDbi")
if(!require("org.Hs.eg.db",quietly=TRUE))
  BiocManager::install("org.Hs.eg.db")
if(!require("biomaRt",quietly=TRUE))
  BiocManager::install("biomaRt")
library(tidyverse)
library(fgsea)
library(AnnotationDbi)
library(org.Hs.eg.db)

getImportantPaths <- function(can){
  print(can)
  allGenes <- read_lines("allGenes.txt")
  print("Got All Genes")
  pathways <- c(gmtPathways("h.all.v2026.1.Hs.entrez.gmt"),gmtPathways("c2.cp.v2026.1.Hs.entrez.gmt"),gmtPathways("c5.all.v2026.1.Hs.entrez.gmt"))
  pathwayNames <- names(pathways)
  allPathways <- pathways
  allPathNames <- pathwayNames
  print("Got Pathways")
  
  allNetGenes <- c()
  allMiRs <- c()
  geneLists <- list()
  miRLists <- list()
  doc <- read.delim(paste(can,"_Modules_Bipartite.txt",sep=""),header=F)
  for(r in rownames(doc)){
    str <- doc[r,]
    vs <- str_split(str,", ")[[1]]
    genes <- vs[!str_detect(vs,"hsa")]
    miRs <- vs[str_detect(vs,"hsa")]
    allNetGenes <- union(allGenes,genes)
    allMiRs <- union(allMiRs,miRs)
    geneLists[[r]] <- genes
    miRLists[[r]] <- miRs
  }
  
  results <- c()
  resultsNoMiR <- list()
  for(r in rownames(doc)){
    print(r)
    modGenes <- as.vector(geneLists[[r]])
    pVals <- c()
    ORs <- c()
    FEs <- c()
    inDF <- c()
    
    pathwayGeneLists <- c()
    
    for(p in allPathNames){
      if(match(p,allPathNames)%%1000 == 0){print(match(p,allPathNames)/length(allPathNames))}
      pathwayGenes <- intersect(allPathways[[p]],allGenes)
      modAndPathwayGenes <- intersect(pathwayGenes,modGenes)
      #print(p)
      #print(modAndPathwayGenes)
      if(length(modAndPathwayGenes) >= 0){
        pathwayGeneLists <- c(pathwayGeneLists,paste(modAndPathwayGenes, collapse=", "))
        a <- length(intersect(modGenes,pathwayGenes))
        b <- length(setdiff(pathwayGenes,modGenes))
        c <- length(setdiff(modGenes,pathwayGenes))
        d <- length(setdiff(allGenes,union(modGenes,pathwayGenes)))
        
        fisher <- fisher.test(rbind(c(a,b),c(c,d)))
        pVals <- c(pVals,fisher$p.value)
        ORs <- c(ORs,fisher$estimate[["odds ratio"]])
        FEs <- c(FEs,(a*length(allGenes)/length(modGenes)/length(pathwayGenes)))
        inDF <- c(inDF,TRUE)
      }
      else{
        inDF <- c(inDF,FALSE)
      }
    }
    padj <- p.adjust(pVals,method="BH")
    fisherDF <- data.frame(allPathNames[inDF],ORs,FEs,pVals,padj,pathwayGeneLists)
    write.csv(fisherDF,paste(can,"_Module",r,"_pathwayAll_FisherResults.csv",sep=""))
    if(nrow(fisherDF) == 0){
      topPaths <- ""
    }
    else{
      topPaths <- (fisherDF %>% filter(FEs > 1, padj < 0.05))
      topPaths <- topPaths[order(topPaths$FEs,decreasing = T),]
      topPaths <- unlist(lapply(rownames(topPaths),function(x){paste(topPaths[x,]$allPathNames.inDF.,"-FE:",topPaths[x,]$FEs,sep = "")}))
    }
    results <- c(results,paste(c(paste("Module",r,sep=""),topPaths),collapse=","))
  }
}
