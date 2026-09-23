# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:06:58 2025

@author: Eli
"""

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import condor
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
    
def getCommunitySignificance(G,modules,plot=False,num_of_rand_net=1000):

    sg, p_values, q_vals, s_vals = qs.qstest(G, modules, qs.qmod, qs.n, condorBRIM, num_of_rand_net=num_of_rand_net, num_of_thread = 1)

    if(plot):
        boxPlotVals = []
        qList = []
        for comNum in range(len(modules)):

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
         
        plt.figure(0,dpi = 200,figsize=(10,5))
        plt.boxplot(boxPlotVals,label="Random Graphs")
        plt.scatter(xSig,qSig,label = "Significant")
        plt.scatter(xNotSig,qNotSig,label = "Not Significant",c='r')
        plt.xlabel("Community",fontsize=12)
        plt.ylabel("Modularity",fontsize=12)
        plt.legend(bbox_to_anchor = (1,1),fontsize=12)
        
        meds = []
        sXs = []
        for s in set(s_vals):
            sXs.append(s)
            meds.append(np.median(q_vals[np.where(s_vals == s)]))
        
        regMed=LinearRegression().fit(np.array(sXs).reshape(-1,1),np.array(meds).reshape(-1,1))
        
        plt.figure(1,dpi=200)
        plt.scatter(sXs,meds,label="Median Modularity of Random Networks")
        plt.scatter([len(m) for m in modules],qList,label="Modularities")
        plt.plot([0,max([len(m) for m in modules])+25],regMed.predict(np.array([0,max([len(m) for m in modules])+25]).reshape(-1,1)),color='k',ls="dashed",label="Linear Regression")
        plt.xlabel("Size",fontsize=12)
        plt.ylabel("Modularity",fontsize=12)
        plt.legend(fontsize=12)
        plt.show()
    return sg, p_values