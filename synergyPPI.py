# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 13:23:45 2026

@author: Eli
"""

import pandas as pd
import networkx as nx
import numpy as np
import mygene
import itertools as itr
from ResistanceCalculation import getResistance
import scipy.stats as sp
mg = mygene.MyGeneInfo()
            
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
            #print(deg)
            for d in range(max(0,deg),deg+1):
                if(d in degree_to_nodes):
                    degNodes = degree_to_nodes[d].copy()
                    if(d == deg):
                        degNodes.remove(node)
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
    
def calculate_closest_distance_median(G, sources, targets, distDict):
    vals = []
    numTargets = 0
    for n1 in targets:
        minDist = np.inf
        distIsFound = False
        for n2 in sources:
            if(n2 in distDict[n1]):
                distIsFound = True
                pathLength = distDict[n1][n2]
                if(pathLength < minDist):
                    minDist = pathLength
        if(distIsFound):
            numTargets += 1
            vals.append(minDist)
    if(numTargets == 0):
        return np.nan
    return np.median(vals)

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
    pathGenes = list(set([x for x in pathList[2:]]))
    if(len(pathGenes) > 0):
        paths.append(pathList[0])
        pathDict[pathList[0]] = pathGenes
        
pathDict["ERK_MAPK_PATHWAY"] = list(set(pathDict["BIOCARTA_ERK_PATHWAY"])|set(pathDict["BIOCARTA_MAPK_PATHWAY"]))


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
PPI = nx.from_pandas_edgelist(edgeDF,"entrez1","entrez2")
"""
cancer = "BRCA"
 
PPI = nx.read_graphml("globalPPI.graphml")

f = open(cancer+"_Genes.txt","r")
lines = f.readlines()
f.close()
topMRNAs = [x.strip() for x in lines]

PPI = nx.subgraph(PPI,topMRNAs)
LCC = sorted(list(nx.connected_components(PPI)),key=len,reverse=True)[0]
PPI = nx.subgraph(PPI, LCC)

print("Got PPI")

resDict = getResistance(PPI)
print("Got Resistance Dict")

G = nx.read_edgelist(cancer+"_EdgeList.txt")
LCC = sorted(list(nx.connected_components(G)),key=len,reverse=True)[0]
G = nx.subgraph(G, LCC)
    
allMRNAs = [n for n in G if "hsa" not in n]
allMiRs = [n for n in G if "hsa" in n]
    
f = open(cancer+"_Modules_Bipartite.txt","r")
modules = f.readlines()
f.close()
modules = [m.rstrip().split(", " ) for m in modules]

indexNames = ["Individual Num Targets",
              "Combined Num Targets",
              "Amount of Overlap",
              "Amount of Pathway Targeted",
              "Individual % Pathway Targeted",
              "Combined % Pathway Targeted",
              "Overlap % Pathway Targeted",
              "Off Target %",
              "Average PPI Resistance From All Targets",
              "PPI Resistance Z-Score (All Targets)",
              "Effectiveness",
              "Effectiveness Z-Score",
              "Effectiveness Test Statistic",
              "Effectiveness p-Value",
              "Precision",
              ]

targetPathways = ["HALLMARK_DNA_REPAIR",
"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
"ERK_MAPK_PATHWAY",
"HALLMARK_PI3K_AKT_MTOR_SIGNALING",
"REACTOME_ESR_MEDIATED_SIGNALING"]


numSamples = 100

for targetedPathway in targetPathways:
    print(targetedPathway)
    resDF = pd.DataFrame(index = indexNames)
    pathGenes = list(set(pathDict[targetedPathway])&set(PPI.nodes()))
    if(len(pathGenes) == 0):
        continue
    for modNum in range(len(modules)):
        print(f"Modules #{modNum}")
        modGenes = [n for n in modules[modNum] if "hsa" not in n]
        modMiRs = [n for n in modules[modNum] if "hsa" in n]
        for miRs in itr.combinations(modMiRs, 2):
            individualMiRTargets = [list(set(list(G.neighbors(x)))&set(PPI.nodes())) for x in miRs]
            nodes_from_random_individual = [get_random_nodes(PPI, x, n_random = numSamples) for x in individualMiRTargets]
            individualTargetNums = [len(x) for x in individualMiRTargets]
            allMiRTargets = list(set([x for x in list(nx.node_boundary(G,miRs))])&set(PPI))
            pctTargeted = len(set(allMiRTargets)&set(pathGenes))/len(pathGenes)
            if(pctTargeted > 0.1):
                overlap = list(set.intersection(*[set(l) for l in individualMiRTargets])&set(PPI))
                numTargeted = len(set(allMiRTargets)&set(pathGenes))
                pctIndividualTargeted = [len(set(x)&set(pathGenes))/len(pathGenes) for x in individualMiRTargets]
                
            
                pctIndividualTargeted_random = [[len(set(x)&set(pathGenes))/len(pathGenes) for x in rand] for rand in nodes_from_random_individual]
                pctOffTarget = 1-(len(set(allMiRTargets)&set(pathGenes))/len(allMiRTargets))
                pctOverlap = len(set(overlap)&set(pathGenes))/len(pathGenes)
                
                dist = calculate_closest_distance_median(PPI, allMiRTargets, pathGenes, resDict)
                
                nodes_from_random_all = get_random_nodes(PPI, allMiRTargets, n_random = numSamples)
                
                nodes_to_random = get_random_nodes(PPI, pathGenes, n_random = numSamples)
                random_values_list = zip(nodes_from_random_all, nodes_to_random)
                vals = np.empty(len(nodes_from_random_all))
                for i, values_random in enumerate(random_values_list):
                    nodes_from, nodes_to = values_random
                    vals[i] = calculate_closest_distance_median(PPI, nodes_from, nodes_to, resDict)
                m, s = np.mean(vals), np.std(vals)
                if(s == 0):
                    z = 0.0
                else:
                    z = (dist-m)/s
                
                    
                z_random = []
                for i in range(len(nodes_from_random_individual[0])):
                    targets = list(set(nodes_from_random_individual[0][i])|set(nodes_from_random_individual[1][i]))
                    d = calculate_closest_distance_median(PPI, targets, pathGenes, resDict)
                    if(s==0):
                        z_random.append(0)
                    else:
                        z_random.append((d-m)/s)
                
                effectiveness = -1*sum(pctIndividualTargeted)*z
                effectivenessRand = [-1*(pctIndividualTargeted_random[0][i]+pctIndividualTargeted_random[1][i])*z_random[i] for i in range(len(z_random))]
                effectivenessZ = (effectiveness-np.mean(effectivenessRand))/np.std(effectivenessRand)
                effectiveness_stat,effectiveness_p = sp.wilcoxon(np.array(effectivenessRand)-effectiveness,alternative="less")
                precision = 1-pctOffTarget
                
                resDF = pd.concat([resDF,pd.DataFrame([str(individualTargetNums)[1:-1],
                                                       len(allMiRTargets),
                                                       len(overlap),
                                                       numTargeted,
                                                       str(pctIndividualTargeted)[1:-1],
                                                       pctTargeted,
                                                       pctOverlap,
                                                       pctOffTarget,
                                                       dist,
                                                       z,
                                                       effectiveness,
                                                       effectivenessZ,
                                                       effectiveness_stat,
                                                       effectiveness_p,
                                                       precision,
                                                       ],
                                                       columns=[", ".join(miRs)],index=indexNames)],axis=1)
    A = resDF.T
    padj = sp.false_discovery_control(list(A["Effectiveness p-Value"].values))
    A.insert(13,"Effectiveness pAdj",padj)
    A = A.sort_values(by='Effectiveness',ascending=False)
    print(targetedPathway)
    print(A[["Effectiveness","Effectiveness pAdj"]].head(3))
    print()
    
