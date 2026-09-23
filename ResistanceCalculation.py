# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:24:42 2026

@author: Eli
"""

import numpy as np
import networkx as nx
import itertools as itr


def getResistance(G,nodes_from=None,nodes_to=None):
    if(nodes_from == None):
        nodes_from = list(G.nodes())
    if(nodes_to == None):
        nodes_to = list(G.nodes())
    resDict = {}
    ccs = list(nx.connected_components(G))
    for i,cc in enumerate(ccs):
        if(len(cc) == 1):
            resDict.update({list(cc)[0]:{list(cc)[0]:0.0}})
            continue
        if(len(cc) == 2):
            nodes = list(cc)
            resDict.update({nodes[0]:{nodes[0]:0.0,nodes[1]:1.0},nodes[1]:{nodes[0]:1.0,nodes[1]:0.0}})
            continue
        subG = nx.subgraph(G,cc)
        nodes = list(subG.nodes())
        nodes_from_cc = list(set(nodes_from)&set(nodes))
        nodes_from_idx = [nodes.index(n) for n in nodes_from_cc]
        nodes_to_cc = list(set(nodes_to)&set(nodes))
        nodes_to_idx = [nodes.index(n) for n in nodes_to_cc]
        
        #A: row normalized adjacency matrix (i.e., random walk matrix)
        A = nx.to_numpy_array(subG)
        rowSums = np.sum(A,axis=0)
        A = np.divide(A,rowSums[:,np.newaxis])
        
        #H: hitting time matrix
        #Need to calculate for each row/column because H(i,i) = 0, so we remove that row/column and calculate all others 
        #H = (I-A)^-1*1Mat
        pct = 1
        C = np.zeros((len(nodes_from_cc),len(nodes_to_cc)))
        for j in range(len(nodes_to_idx)):#range(len(A)):
            if(j/(len(nodes_from_idx)+len(nodes_to_idx))*100>pct):
                print(pct)
                pct += 1
            P = np.delete(A,nodes_to_idx[j],axis=0)
            P = np.delete(P,nodes_to_idx[j],axis=1)
            I = np.identity(len(P))
            inv = np.linalg.inv(I-P)
            eta = inv@np.ones((len(P),1))
            etaRes = np.insert(eta,nodes_to_idx[j],0)
            C[:,j] += etaRes[nodes_from_idx]
        for j in range(len(nodes_from_idx)):#range(len(A)):
            if((j+len(nodes_to_idx))/(len(nodes_from_idx)+len(nodes_to_idx))*100>pct):
                print(pct)
                pct += 1
            P = np.delete(A,nodes_from_idx[j],axis=0)
            P = np.delete(P,nodes_from_idx[j],axis=1)
            I = np.identity(len(P))
            inv = np.linalg.inv(I-P)
            eta = inv@np.ones((len(P),1))
            etaRes = np.insert(eta,nodes_from_idx[j],0)
            C[j,:] += etaRes[nodes_to_idx]
        
        #print(H)
        #C: Commute Time = H(i,j)+H(j,i)
        #C = H+H.T
        #print(C)
        #R: Resistance, C = 2*m*Resistance
        R = C/2/len(subG.edges())
        #print(R)
        subDict = {nodes_from_cc[i]:{nodes_to_cc[j]:R[i,j] for j in range(len(nodes_to_idx))} for i in range(len(nodes_from_idx))}
        resDict.update(subDict)
    return resDict