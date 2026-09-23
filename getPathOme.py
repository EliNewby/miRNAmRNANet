# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:57:47 2026

@author: Eli
"""

import networkx as nx

def getPathNet(paths,pathDict,jaccardThreshold=0.0,verbose=True):
    pathNet = nx.Graph()
    pathNet.add_nodes_from(paths)
    pct = 1
    edges = []
    overlaps = []
    unions = []
    for i in range(len(paths)):
        p = paths[i]
        if(verbose):
            if(i/len(paths)*100 > pct):
                print(pct)
                pct += 1
        for j in range(i):
            p2 = paths[j]
            edges.append((p,p2))
            overlap = len(set(pathDict[p])&set(pathDict[p2]))
            overlaps.append(overlap)
            unions.append(len(set(pathDict[p])|set(pathDict[p2])))
    #pAdjs = sp.false_discovery_control(pVals)
    jaccards = [overlaps[i]/unions[i] for i in range(len(overlaps))]
    
    edgesSig = [edges[i]+(jaccards[i],) for i in range(len(edges)) if jaccards[i] > jaccardThreshold]
    pathNet.add_weighted_edges_from(edgesSig)
    return pathNet
