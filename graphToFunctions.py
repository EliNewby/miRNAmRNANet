# -*- coding: utf-8 -*-
"""
Created on Tue May 20 14:02:01 2025

@author: Eli
"""

import rpy2
import os
os.environ['R_HOME'] = "C:/Program Files/R/R-4.3.3"
import rpy2.robjects as ro
import networkx as nx
import pandas as pd
import numpy as np
import itertools as itr
import condor
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import random
import scipy.stats as sp
from validateCommunities import getCommunitySignificance
from ResistanceCalculation import getResistance
from getPathOme import getPathNet
import sys
#import mygene
#mg = mygene.MyGeneInfo()
"""
df = pd.read_csv("C:/Users/Eli/Documents/Research/CleClinic/9606.protein.info.v12.0.txt",sep="\t")
nameDict = {df.loc[i,'#string_protein_id']:df.loc[i,'preferred_name'] for i in range(len(df))}

A=mg.querymany(list(nameDict.values()),scopes="symbol",fields="entrezgene",species="human",returnall=True)
miss = A["missing"]
dups = [a[0] for a in A["dup"]]

symbolToEntrez = {}
for a in A["out"]:
    if(a["query"] in miss):
        continue
    if(a["query"] in dups):
        if("entrezgene" in a):
            symbolToEntrez.update({a["query"]:a["entrezgene"]})
    else:
        symbolToEntrez.update({a["query"]:a["entrezgene"]})

entrezDict = {k:symbolToEntrez[v] for k,v in nameDict.items() if v in symbolToEntrez}

edgeDF = pd.read_csv("C:/Users/Eli/Documents/Research/CleClinic/9606.protein.links.v12.0.onlyAB.txt",sep=" ")
edgeDF = edgeDF[edgeDF["combined_score"]>=700]
protein1 = list(edgeDF["protein1"])
protein2 = list(edgeDF["protein2"])
edgeDF = edgeDF[[protein1[i] in entrezDict and protein2[i] in entrezDict for i in range(len(edgeDF))]]
edgeDF["entrez1"] = [entrezDict[x] for x in edgeDF["protein1"]]
edgeDF["entrez2"] = [entrezDict[x] for x in edgeDF["protein2"]]
globalPPI = nx.from_pandas_edgelist(edgeDF,"entrez1","entrez2")
"""
globalPPI = nx.read_graphml("globalPPI.graphml")

random.seed(5151)
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors2 = list(mcolors.XKCD_COLORS.values())
random.shuffle(colors2)
colors += colors2

rSource = ro.r('''
               source('analyzeModulesFunction.R')
               ''')
r_getImportantPaths = ro.globalenv["getImportantPaths"]

def getModules(cancerType,numCommItrs=100,communityInclusionThreshold=0.75):
    G = nx.read_edgelist(cancerType+"_EdgeList.txt")
    #G = nx.read_edgelist("C:/Users/Eli/Downloads/BRCA_EdgeList_matrix.txt")
    LCC = sorted(list(nx.connected_components(G)),key=len,reverse=True)[0]
    G = nx.subgraph(G, LCC)
    mRNAs = []
    miRs = []
    for n in G:
        if("hsa" in n):
            miRs.append(n)
        else:
            mRNAs.append(n)
    edgeList = []
    for miR in miRs:
        for mRNA in G.neighbors(miR):
            edgeList.append([miR,mRNA])
    edgeList = pd.DataFrame(edgeList)
    import random
    commNet = nx.Graph()
    for i in range(numCommItrs):
        random.seed(i)
        co = condor.condor_object(dataframe=edgeList,silent=True)
        co.initial_community()
        co.brim()
        
        modDict = {}
        for i in range(len(co.tar_memb)):
            if(co.tar_memb.iloc[i]["community"] not in modDict):
                modDict[co.tar_memb.iloc[i]["community"]] = [co.tar_memb.iloc[i]["tar"][4:]]
            else:
                modDict[co.tar_memb.iloc[i]["community"]].append(co.tar_memb.iloc[i]["tar"][4:])
        for i in range(len(co.reg_memb)):
            if(co.reg_memb.iloc[i]["community"] not in modDict):
                modDict[co.reg_memb.iloc[i]["community"]] = [co.reg_memb.iloc[i]["reg"][4:]]
            else:
                modDict[co.reg_memb.iloc[i]["community"]].append(co.reg_memb.iloc[i]["reg"][4:])
        modDict = {k:modDict[k] for k in sorted(modDict)}
        for k in modDict:
            comm = modDict[k]
            for i1 in range(len(comm)):
                n1 = comm[i1]
                for i2 in range(i1):
                    n2 = comm[i2]
                    if(commNet.has_edge(n1,n2)):
                        commNet.edges()[(n1,n2)]["numItrs"] += 1
                    else:
                        commNet.add_edge(n1,n2,numItrs = 1)
                        
    consensusEdges = [x for x in commNet.edges() if commNet.edges()[x]["numItrs"] > communityInclusionThreshold*numCommItrs]
    consensusNet = commNet.edge_subgraph(consensusEdges)
    ccs = list(nx.connected_components(consensusNet))
    modules = []
    for cc in ccs:
        hasMiR = len([x for x in cc if "hsa" in x]) > 0
        hasMRNA = len([x for x in cc if "hsa" not in x]) > 0
        if(hasMiR and hasMRNA):
            modules.append(cc)
    modules = [sorted(x) for x in modules if len(x) > 0]
    sig,pVals = getCommunitySignificance(G, modules,plot=True)

    f = open(f"{cancerType}_Modules_Bipartite.txt",'w')
    for i in range(len(modules)):
        if(sig[i]):
            string = ", ".join(modules[i])
            string += "\n"
            f.write(string)
    f.close()

    
def get_random_nodes(G, nodes_selected, n_random, seed=None):
    if(seed is not None):
        rng = np.random.default_rng(seed)
    else:
        rng = np.random.default_rng()
    degree_to_nodes = {}
    for node, degree in G.degree():
        degree_to_nodes.setdefault(degree, []).append(node)

    randNodes = []
    for i in range(n_random):
        nodes_random = []
        for node in nodes_selected:
            deg = G.degree(node)
            equivalentNodes = []
            for d in range(max(0,deg),deg+1):
                if(d in degree_to_nodes):
                    degNodes = degree_to_nodes[d].copy()
                    equivalentNodes += degNodes
            if(len(equivalentNodes) < 10):
                allDegs = np.array(list(degree_to_nodes.keys()))
                distFromDeg = abs(allDegs-deg)
                closestDegs = sorted(distFromDeg)[1:]
                i = 0
                while(len(equivalentNodes) < 10):
                    dist = closestDegs[i]
                    i += 1
                    closeDeg = allDegs[np.where(distFromDeg == dist)[0][0]]
                    equivalentNodes += degree_to_nodes[closeDeg]
            chosen = rng.choice(equivalentNodes)
            for k in range(20):
                if(chosen in nodes_random):
                    chosen = rng.choice(equivalentNodes)
            nodes_random.append(chosen)
        randNodes.append(nodes_random)
    return randNodes

def calculate_closest_distance(G, sources, targets, distDict):
    minSum = 0
    numTargets = 0
    for n1 in targets:
        minDist = 10000000
        for n2 in sources:
            pathLength = distDict[n1][n2]
            if(pathLength < minDist):
                minDist = pathLength
    return minSum/len(targets)

def calculate_closest_distance_median(G, sources, targets, distDict):
    vals = []
    numTargets = 0
    for n1 in targets:
        minDist = 10000000
        for n2 in sources:
            pathLength = distDict[n1][n2]
            if(pathLength < minDist):
                minDist = pathLength
        vals.append(minDist)
    return np.median(vals)

def getPPIProximityFunctions(cancerType, numSamples = 1000):
    G = nx.read_edgelist(cancerType+"_EdgeList.txt")
    LCC = sorted(list(nx.connected_components(G)),key=len,reverse=True)[0]
    G = nx.subgraph(G, LCC)
    
    f = open(cancerType+"_Genes.txt","r")
    lines = f.readlines()
    f.close()
    topMRNAs = [x.strip() for x in lines]
    
    
    PPI = nx.subgraph(globalPPI,topMRNAs)
    LCC = sorted(list(nx.connected_components(PPI)),key=len,reverse=True)[0]
    PPI = nx.subgraph(PPI, LCC)
    
    resisDict = getResistance(PPI)
      
    print("Got Distance Dict")
    
    f = open(cancerType+"_Modules_Bipartite.txt","r")
    modules = f.readlines()
    f.close()
    modules = [m.rstrip().split(", " ) for m in modules]
    significantPaths = set()
    modRandomNodesDict = {}
    for modNum in range(len(modules))[7:]:
        print(f"Module {modNum}")
        modPathDF =  pd.read_csv(cancerType+"_Module"+str(modNum+1)+"_pathwayAll_FisherResults.csv",index_col=0)
        sigPathDF = modPathDF[modPathDF["padj"] < 0.05]
        significantPaths |= set(sigPathDF["allPathNames.inDF."])
        modMiRs = [n for n in modules[modNum] if "hsa" in n]
        miRTargets = [x for x in list(nx.node_boundary(G,modMiRs))]
        miRTargets = set(miRTargets) & set(PPI.nodes())
        modRandomNodesDict[modNum] = get_random_nodes(PPI, miRTargets, n_random = numSamples)
    significantPaths = list(significantPaths)
    #significantPaths = ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
    print("Got Significant Paths + miR Target Random Nodes")
    lines = []
    f = open("h.all.v2026.1.Hs.entrez.gmt","r")
    lines += f.readlines()
    f.close()
    
    f = open("c2.cp.v2026.1.Hs.entrez.gmt","r")
    lines += f.readlines()
    f.close()
    
    f = open("c5.all.v2026.1.Hs.entrez.gmt","r")
    lines += f.readlines()
    f.close()
    
    paths = []
    pathDict = {}
    pathRandomNodesDict = {}
    pct = 1
    for i in range(len(lines)):
        if(i/len(lines)*100 > pct):
            print(f"{pct}%")
            pct += 1
        l = lines[i].rstrip()
        pathList = l.split('\t')
        pathGenes = list(set([x for x in pathList[2:]])&set(PPI.nodes()))
        if(len(pathGenes) > 0 and pathList[0] in significantPaths):
            paths.append(pathList[0])
            pathDict[pathList[0]] = pathGenes
            pathRandomNodesDict[pathList[0]] = get_random_nodes(PPI, pathGenes, n_random = numSamples)
    print("Got Paths + Path Random Nodes")
    for modNum in list(range(len(modules)))[7:]:
        print(f"Modules {modNum}")
        modMiRs = [n for n in modules[modNum] if "hsa" in n]
        miRTargets = [x for x in list(nx.node_boundary(G,modMiRs))]
        miRTargets = set(miRTargets) & set(PPI.nodes())
        
        nodes_from_random = modRandomNodesDict[modNum]
        modPathDF =  pd.read_csv(cancerType+"_Module"+str(modNum+1)+"_pathwayAll_FisherResults.csv",index_col=0)
        sigPathDF = modPathDF[modPathDF["padj"] < 0.05]
        modSigPaths = list(sigPathDF["allPathNames.inDF."])
        #modSigPaths = ['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']  
        pathZs = []
        pathZs_pVals = []
        pathZMods = []
        dists = []
        medians = []
        MADs = []
        pct = 1
        sigPaths = []
        for pathNum,p in enumerate(modSigPaths):
            if(pathNum/len(modSigPaths)*100 > pct):
                print(str((pathNum*100)//len(modSigPaths))+"%")
                pct = (pathNum*100)//len(modSigPaths) + 1
            if(p not in pathDict):
                continue
            sigPaths.append(p)
            pathGenes = pathDict[p]
            #d = calculate_closest_distance(PPI, miRTargets, pathGenes, distDict)
            d = calculate_closest_distance_median(PPI, miRTargets, pathGenes, resisDict)
            dists.append(d)
            
            nodes_to_random = pathRandomNodesDict[p]
            random_values_list = zip(nodes_from_random, nodes_to_random)
            values = np.empty(len(nodes_from_random))
            for i, values_random in enumerate(random_values_list):
                nodes_from, nodes_to = values_random
                #values[i] = calculate_closest_distance(PPI, nodes_from, nodes_to, distDict)
                values[i] = calculate_closest_distance_median(PPI, nodes_from, nodes_to, resisDict)
            
            mean = np.mean(values)
            std = np.std(values)
            if(std == 0):
                z = 0.0
            else:
                z = (d-mean)/std
            pathZs.append(z)    
            pathZs_pVals.append(1-sp.norm.sf(z))
            
            med = np.median(values)
            medians.append(med)
            MAD = np.median(abs(values-med))
            MADs.append(MAD)
            if(MAD == 0):
                zMod = 0.0
            else:
                zMod = 0.6745*(d-med)/MAD
            
            pathZMods.append(zMod)
        print("Got Significant Pathway, PPI Z-Scores")
        resDF = pd.DataFrame([pathZs,pathZs_pVals,pathZMods,dists,medians,MADs],columns=sigPaths,index=["Z-Score","pVals","Modified Z-Score","Distance","Median","MAD"]).transpose()
        pAdj = sp.false_discovery_control(resDF["pVals"])
        resDF.insert(2,column="pAdj",value=pAdj)
        resDF = resDF.sort_values(by="Modified Z-Score")
        resDF.to_excel(cancerType+"_pathwayZScores_Module"+str(modNum)+".xlsx")

def getModuleFunctions(cancerType):
    G = nx.read_edgelist(cancerType+"_EdgeList.txt")
    LCC = sorted(list(nx.connected_components(G)),key=len,reverse=True)[0]
    G = nx.subgraph(G, LCC)
    
    lines = []
    f = open("h.all.v2026.1.Hs.entrez.gmt","r")
    lines += f.readlines()
    f.close()
    
    f = open("c2.cp.v2026.1.Hs.entrez.gmt","r")
    lines += f.readlines()
    f.close()
    
    f = open("c5.all.v2026.1.Hs.entrez.gmt","r")
    lines += f.readlines()
    f.close()
    
    paths = []
    pathDict = {}
    for i in range(len(lines)):
        l = lines[i].rstrip()
        pathList = l.split('\t')
        pathGenes = [x for x in pathList[2:]]
        if(len(pathGenes) > 0):
            paths.append(pathList[0])
            pathDict[pathList[0]] = pathGenes
    
    f = open(cancerType+"_Modules_Bipartite.txt","r")
    modules = f.readlines()
    f.close()
    modules = [m.rstrip().split(", " ) for m in modules]
    writer = pd.ExcelWriter(cancerType+"_moduleFunctions.xlsx",engine="xlsxwriter")
    for modNum in range(len(modules)):
        print(modNum)
        
        targetPathways =  pd.read_csv(cancerType+"_Module"+str(modNum+1)+"_pathwayAll_FisherResults.csv",index_col=0)
        if(len(targetPathways)==0):
            continue
        relevantTargetPaths = list(targetPathways[(targetPathways["FEs"] > 0) & (targetPathways["padj"] < 0.05)]["allPathNames.inDF."])
        
        protDF = pd.read_excel(cancerType+"_pathwayZScores_Module"+str(modNum)+".xlsx")
        relevantProtPaths = list(protDF[(protDF["pAdj"] < 0.05)]["Unnamed: 0"])
        
        relevantPaths = list(set(relevantTargetPaths)&set(relevantProtPaths))
        relevantPaths = [x for x in relevantPaths if x in paths]
        uniquePaths = set(relevantPaths)
        
        for i in range(len(modules)):
            if(i != modNum):
                targetPathways2 =  pd.read_csv(cancerType+"_Module"+str(i+1)+"_pathwayAll_FisherResults.csv",index_col=0)
                if(len(targetPathways2) == 0):
                    continue
                relevantTargetPaths2 = list(targetPathways2[(targetPathways2["FEs"] > 0) & (targetPathways2["padj"] < 0.05)]["allPathNames.inDF."])
                
                protDF2 = pd.read_excel(cancerType+"_pathwayZScores_Module"+str(i)+".xlsx")
                relevantProtPaths2 = list(protDF2[(protDF2["pAdj"] < 0.05)]["Unnamed: 0"])
                
                relevantPaths2 = list(set(relevantTargetPaths2)&set(relevantProtPaths2))
                uniquePaths = uniquePaths - set(relevantPaths2)
                 
        uniquePaths = list(uniquePaths)
        
        
        pathRanks = {}
        for i in range(len(relevantPaths)):
            p = relevantPaths[i]
            FE = targetPathways[targetPathways["allPathNames.inDF."]==p]["FEs"].values[0]
            prot = protDF[protDF["Unnamed: 0"] == p]["Z-Score"].values[0]
            effectiveness = -1*FE*prot
            pathRanks[p] = effectiveness
        pathRanks = {k:v for k,v in sorted(pathRanks.items(),key=lambda x:x[1],reverse=True)}

        subG = nx.subgraph(G,modules[modNum])
        miRDegs = {n:subG.degree(n) for n in subG if "hsa" in n}
        miRDegs = {k:v for k,v in sorted(miRDegs.items(),key=lambda x:x[1])}
        print("\n".join([f"{k}:{v}" for k,v in miRDegs.items()]))
        print()
        
        pathG = getPathNet(relevantPaths,pathDict,verbose=False,jaccardThreshold=0.1)
        
        for p in pathG:
            pathG.nodes()[p]["effectiveness"] = pathRanks[p]
            
        ccs = nx.community.greedy_modularity_communities(pathG,weight="weight")
        sortedCCs = sorted(ccs,key=lambda x:np.average([pathG.nodes()[n]["effectiveness"] for n in x]),reverse=True)
        resDF = pd.DataFrame(columns = ["FE","FE-pVal","FE-pAdj","Z-Score","Z-pVal","Z-pAdj","isUnique","Pathway Community"])
        for i,cc in enumerate(sortedCCs):
            sortedFuncs = sorted(cc,key=lambda x:pathG.nodes()[x]["effectiveness"],reverse=True)
            for p in sortedFuncs:
                FE = targetPathways[targetPathways["allPathNames.inDF."]==p]["FEs"].values[0]
                FEpVal = targetPathways[targetPathways["allPathNames.inDF."]==p]["pVals"].values[0]
                FEpAdj = targetPathways[targetPathways["allPathNames.inDF."]==p]["padj"].values[0]
                prot = protDF[protDF["Unnamed: 0"] == p]["Z-Score"].values[0]
                protPVal = protDF[protDF["Unnamed: 0"] == p]["p"].values[0]
                protPAdj = protDF[protDF["Unnamed: 0"] == p]["pAdj"].values[0]
                isUnique = p in uniquePaths
                resDF.loc[p] = [FE,FEpVal,FEpAdj,prot,protPVal,protPAdj,isUnique,i]
        resDF.to_excel(writer,sheet_name=f"Module {modNum}")
    writer.close()

cancers = ["BRCA","BRCA_Basal", "BRCA_HER2", "BRCA_LumA", "BRCA_LumB"]
netDict = {i:cancers[i] for i in range(len(cancers))}
for i in range(len(cancers)):#[int(sys.argv[1])]:#cancers:
    cancerType = netDict[i]
    print(cancerType)
    #getModules(cancerType)
    #r_getImportantPaths(cancerType)
    #getPPIProximityFunctions(cancerType,numSamples=1000)
    #getModuleFunctions(cancerType)
    