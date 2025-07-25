# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:06:58 2025

@author: Eli
"""

import networkx as nx
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import itertools as itr
import condor
import matplotlib.colors as mcolors
import random
import qstest as qs
from sklearn.linear_model import LinearRegression

def condorBRIM(network):
    set1,set2 = nx.bipartite.sets(network)
    edgeList = []
    for n1 in set1:
        for n2 in network.neighbors(n1):
            edgeList.append([n1,n2])
    edgeList = pd.DataFrame(edgeList)
    co = condor.condor_object(dataframe=edgeList,silent=True)
    co.initial_community()
    co.brim()
    
    modDict = {}
    for i in range(len(co.tar_memb)):
        if(co.tar_memb.iloc[i]["community"] not in modDict):
            modDict[co.tar_memb.iloc[i]["community"]] = [int(co.tar_memb.iloc[i]["tar"][4:])]
        else:
            modDict[co.tar_memb.iloc[i]["community"]].append(int(co.tar_memb.iloc[i]["tar"][4:]))
    for i in range(len(co.reg_memb)):
        if(co.reg_memb.iloc[i]["community"] not in modDict):
            modDict[co.reg_memb.iloc[i]["community"]] = [int(co.reg_memb.iloc[i]["reg"][4:])]
        else:
            modDict[co.reg_memb.iloc[i]["community"]].append(int(co.reg_memb.iloc[i]["reg"][4:]))
    modDict = {k:modDict[k] for k in sorted(modDict)}
    
    communities = [x for x in modDict.values()]
    return communities
    

#cancers = ["BRCA","LUAD","COAD","OV","BLCA","LGG","KIRC","LUSC","ACC","CESC","ESCA","HNSC","KICH","KIRP","LIHC","MESO","PAAD","PCPG","LAML","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
cancers = ["BRCA"]
#cancers = ["Basal","Her2","LumA","LumB"]
numNets = len(cancers)
netDict = {i:cancers[i] for i in range(len(cancers))}

maxModNum = 0
for netNum in range(numNets):
    f = open(netDict[netNum]+"/"+netDict[netNum]+"_Modules_Bipartite_nearMin.txt","r")
    #f = open("BRCA/BRCA_Modules_Bipartite_"+netDict[netNum]+".txt","r")
    modules = f.readlines()
    f.close()
    if(len(modules) > maxModNum):
        maxModNum = len(modules)

sigDF = pd.DataFrame(columns = cancers, index=np.arange(maxModNum))
pValDF = pd.DataFrame(columns = cancers, index=np.arange(maxModNum))

for netNum in range(numNets):
    print(netDict[netNum])
    #G = nx.read_graphml(netDict[netNum]+'/'+netDict[netNum]+".graphml")
    G = nx.read_edgelist(netDict[netNum]+'/'+netDict[netNum]+"_EdgeList_nearMin.txt")
    #G = nx.read_edgelist("BRCA/BRCA_EdgeList_"+netDict[netNum]+".txt")
    LCC = sorted(list(nx.connected_components(G)),key=len,reverse=True)[0]
    G = nx.subgraph(G, LCC)
    
    f = open(netDict[netNum]+"/"+netDict[netNum]+"_Modules_Bipartite_nearMin.txt","r")
    #f = open("BRCA/BRCA_Modules_Bipartite_"+netDict[netNum]+".txt","r")
    modules = f.readlines()
    f.close()
    modules = [m.rstrip().split(", " ) for m in modules]

    sg, p_values, q_vals, s_vals = qs.qstest(G, modules, qs.qmod, qs.n, condorBRIM, num_of_rand_net=1000, num_of_thread = 1)
    sigDF[netDict[netNum]] = sg+[np.nan]*(maxModNum-len(sg))
    pValDF[netDict[netNum]] = p_values+[np.nan]*(maxModNum-len(p_values))
    
    boxPlotVals = []
    qList = []
    for comNum in range(len(modules)):
        """
        degSum = 0
        for n in modules[comNum]:
            degSum += G.degree(n)
        print(degSum)
        
        boxPlotVals.append(q_vals[np.where(abs(np.array(s_vals)-degSum) < 10)])
        """
        
        boxPlotVals.append(q_vals[np.where(abs(np.array(s_vals)-len(modules[comNum])) < 10)])
        community = modules[comNum]
        deg = G.degree(community)
        q = 0
        D = 0
        for i in community:
            for j in community:
                if G.has_edge(i, j) == False:
                    continue
                q += 1.0
            D += deg[i]

        M = G.size() / 2
        q = (q - D * D / (2.0 * M)) / (2 * M)
        qList.append(q)
    
    """
    y = [np.percentile(b,50) for b in boxPlotVals if len(b) > 0]
    x = [len(modules[i]) for i in range(len(modules)) if len(boxPlotVals[i]) > 0]
    regMed=LinearRegression().fit(np.array(x).reshape(-1,1),np.array(y).reshape(-1,1))
    med = regMed.predict(np.array(len(modules[i])).reshape(-1,1))[0,0]
    """
    """
    print("T-Test")
    for i in range(len(modules)):
        print((i+1),sp.stats.ttest_1samp(np.array(s_vals),np.array(len(modules[i])))[1])
    """
    
    xSig = []
    qSig = []
    xNotSig = []
    qNotSig = []
    for i in range(len(sg)):
        if(sg[i]):
            qSig.append(qList[i])
            xSig.append(i+1)
        else:
            qNotSig.append(qList[i])
            xNotSig.append(i+1)
     
    plt.figure(0,dpi = 200)
    plt.boxplot(boxPlotVals,label="Random Graphs")
    plt.scatter(xSig,qSig,label = "Significant")
    plt.scatter(xNotSig,qNotSig,label = "Not Significant",c='r')
    plt.xlabel("Community")
    plt.ylabel("Modularity")
    plt.legend(bbox_to_anchor = (1,1))
    
    meds = []
    sXs = []
    for s in set(s_vals):
        sXs.append(s)
        meds.append(np.median(q_vals[np.where(s_vals == s)]))
    
    regMed=LinearRegression().fit(np.array(sXs).reshape(-1,1),np.array(meds).reshape(-1,1))
    
    plt.figure(1,dpi=200)
    plt.scatter(sXs,meds,label="Median Modularity of Random Networks")
    plt.scatter([len(m) for m in modules],qList,label="BRCA Modularities")
    plt.plot([0,max([len(m) for m in modules])+25],regMed.predict(np.array([0,max([len(m) for m in modules])+25]).reshape(-1,1)),color='k',ls="dashed",label="Linear Regression")
    plt.xlabel("Size")
    plt.ylabel("Modularity")
    plt.legend()
    plt.show()
    
#writer = pd.ExcelWriter('BRCA/SubtypeModuleSignificance2.xlsx', engine='xlsxwriter')
#sigDF.to_excel(writer,sheet_name="Significance")
#pValDF.to_excel(writer,sheet_name="p-values")
#writer.close()