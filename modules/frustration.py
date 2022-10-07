# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 12:52:10 2020

@author: Namhee Kim
"""
import copy
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
import itertools
from pyboolnet.interaction_graphs import primes2igraph
from pyboolnet.file_exchange import primes2bnet
from pyboolnet.prime_implicants import create_constants, percolate, find_constants


# ================================================================================
# compute 'molecular state ambiguity (delta r)'
# by calculating the ratio of frustrated_edge from perturbed(new) truth table
# ================================================================================

# Function for splitting graph
def splitG(primes, fixed_nodes):
    graph = primes2igraph(primes)
    fix_subgraph = graph.subgraph(fixed_nodes)
    nodeList = list(primes.keys())
    
    Adj = nx.adjacency_matrix(graph, nodelist = list(primes.keys())) #source X target
    res_node_idx = [nodeList.index(x) for x in set(nodeList)-fixed_nodes]
    res_idx_out = set(itertools.chain(*[np.where(x)[0] for x in (Adj[res_node_idx,:].toarray() != 0)]))
    res_idx_in = set(itertools.chain(*[np.where(x)[0] for x in (Adj[:,res_node_idx].toarray().T != 0)]))
    res_idx = set(res_node_idx) | res_idx_in | res_idx_out
    
    # The residual_graph consists of edges which are not affected by perturbation
    residual_subgraph = graph.subgraph([nodeList[x] for x in res_idx])
    res_subgraph = nx.DiGraph()
    res_subgraph.add_edges_from([x for x in residual_subgraph.edges(data=True) if x not in list(fix_subgraph.edges(data=True))])
    
    return Adj, graph, fix_subgraph, res_subgraph

def get_bin(x, xlen):
    gb = format(x,'b')
    return ''.join(['0' for x in range(xlen-len(gb))]) + gb


# Function for reducing primes
def reduce_primes(primes: dict, constants: dict):
    reduced_primes = create_constants(primes, constants, copy=True)
    logic_reduced_primes = percolate(primes=reduced_primes, remove_constants=False, copy=True, max_iterations=len(primes)*2)
    fix_dict_reduced = find_constants(logic_reduced_primes)

    logic_reduced_primes = percolate(primes=reduced_primes, remove_constants=True, copy=True, max_iterations=len(primes)*2)  
    if logic_reduced_primes == {}: logic_reduced = ''
    else: logic_reduced = sorted([x for x in primes2bnet(logic_reduced_primes).split('\n') if x != ''])
    
    return logic_reduced_primes, fix_dict_reduced, logic_reduced



def computeFusingT_updated(allperturbs, fix_dict, ctrl, primes, model_file, saveFilename):
    frustrated_out = pd.DataFrame([])
    frustrated_in = pd.DataFrame([])
    allperturbs.append(ctrl)
    allperturbs = list(map(dict, set(tuple(sorted(d.items())) for d in allperturbs)))
    for perturbs in allperturbs:

        # make residual_graph and fixed_graph to categorize edge conditions from the truth table
        fix_dict_tmp = copy.deepcopy(fix_dict)
        fix_dict_tmp.update(perturbs)
        reduced_prime, fix_dict_reduced, _ = reduce_primes(primes, fix_dict_tmp)
        fix_dict_reduced = {x:[-1 if y==0 else 1][0] for x,y in fix_dict_reduced.items()} # for computing frustrated edge
        fixed_nodes = set(fix_dict_reduced.keys())
        Adj, graph, fixed_graph, residual_graph = splitG(primes, fixed_nodes) # The residual_graph consists of edges which are not affected by perturbation
        
        resN = set(reduced_prime.keys()) 
        res_fixN = set(residual_graph.nodes)- set(reduced_prime.keys())
        fixN = set(fixed_graph.nodes) - res_fixN #only in the fixed_graph
        nodeList = list(primes.keys())
        
        # make new trtuh table using reduced_prime after perturbation(s) 
        def reducedprime2T(Adj, reduced_prime, node): 
            allnodes = [nodeList[x] for x in np.where(Adj[:,nodeList.index(node)].toarray().T[0]!=0)[0]]
            fixed = {x:[0 if y==-1 else y][0] for x,y in fix_dict_reduced.items() if x in allnodes}                
            prime0 = reduced_prime[node][0]    
            newprime0 = []
            for p0 in prime0:
                fixedtmp = copy.deepcopy(fixed)
                fixedtmp.update(p0)       
                wcname = list(set(allnodes)-set(fixedtmp.keys()))
                if len(wcname)>0:
                    states = [{x:int(y) for x,y in zip(wcname, get_bin(w, len(wcname)))} for w in range(2**len(wcname))]
                    for sts in states:
                        fixedtmps = copy.deepcopy(fixedtmp)
                        fixedtmps.update(sts)
                        newprime0.append(dict(sorted(fixedtmps.items())))                
                else: newprime0.append(dict(sorted(fixedtmp.items())))
                    
            prime1 = reduced_prime[node][1]  
            newprime1 = []
            for p1 in prime1:
                fixedtmp = copy.deepcopy(fixed)
                fixedtmp.update(p1)       
                wcname = list(set(allnodes)-set(fixedtmp.keys()))
                if len(wcname)>0:
                    states = [{x:int(y) for x,y in zip(wcname, get_bin(w, len(wcname)))} for w in range(2**len(wcname))]
                    for sts in states:
                        fixedtmps = copy.deepcopy(fixedtmp)
                        fixedtmps.update(sts)
                        newprime1.append(dict(sorted(fixedtmps.items())))                
                else: newprime1.append(dict(sorted(fixedtmp.items())))
                    
            newprime0 = list(map(dict, set(tuple(sorted(d.items())) for d in newprime0)))
            newprime1 = list(map(dict, set(tuple(sorted(d.items())) for d in newprime1)))
    
            return {node:[newprime0, newprime1]}
        

        # compute frustrated_edge_ratio of each node
        def computeF_inout(node):
            allF_dict = defaultdict(list)
            if node in resN:
                newt = reducedprime2T(Adj, reduced_prime, node)
                lut_len = len(newt[node][0]) + len(newt[node][1])
                for l0 in newt[node][0]: # false(0)
                    for k,v in l0.items(): allF_dict[k] += [-(1/lut_len)*list(graph[k][node]['sign'])[0]*(-1)*(-1) if v==0 else -(1/lut_len)*list(graph[k][node]['sign'])[0]*1*(-1)] 
                
                for l1 in newt[node][1]: # true(1)
                    for k,v in l1.items(): allF_dict[k] += [-(1/lut_len)*list(graph[k][node]['sign'])[0]*(1)*(-1) if v==0 else (1/lut_len)*list(graph[k][node]['sign'])[0]*1*(-1)] 

            elif node in fixN:
                for source in [nodeList[x] for x in np.where((Adj[:,nodeList.index(node)].toarray().T[0] != 0))[0]]:
                    allF_dict[source] += [fix_dict_reduced[source]*fix_dict_reduced[node]*list(graph[source][node]['sign'])[0]*(-1)]
  
            else:
                for source in [nodeList[x] for x in np.where((Adj[:,nodeList.index(node)].toarray().T[0] != 0))[0]]:
                    if source not in fix_dict_reduced.keys(): 
                        newt = reducedprime2T(Adj, reduced_prime, source)
                        lut_len = len(newt[source][0]) + len(newt[source][1])                        
                        r0 = (-1)*(-1)*fix_dict_reduced[node]*list(graph[source][node]['sign'])[0] * len(newt[source][0])/lut_len
                        r1 = (-1)*(1)*fix_dict_reduced[node]*list(graph[source][node]['sign'])[0] * len(newt[source][1])/lut_len                        
                        allF_dict[source] += [r0+r1]                  
                    else:  allF_dict[source] += [fix_dict_reduced[source]*fix_dict_reduced[node]*list(graph[source][node]['sign'])[0]*(-1)]
                    
            return allF_dict

        # compute frustrated_edge_ratio
        indvTable = pd.DataFrame()
        for node in nodeList:
            tv = computeF_inout(node)
            tv = pd.DataFrame.from_dict(tv, orient = 'index').sum(axis = 1)
            indvTable = pd.concat([indvTable, pd.DataFrame(tv,columns = [node])], axis=1, sort=True)
        indvTable = indvTable.fillna(0)

        f_outinfluence = pd.DataFrame(np.sum(indvTable,axis=1), columns = [str(sorted(perturbs.items()))]) 
        f_ininfluence =  pd.DataFrame(np.sum(indvTable,axis=0), columns = [str(sorted(perturbs.items()))]) 
        frustrated_out = pd.concat([frustrated_out, f_outinfluence], axis=1, sort=True)
        frustrated_in = pd.concat([frustrated_in, f_ininfluence], axis=1, sort=True)

    frustrated_out.loc['Fvalue',:] = np.sum(frustrated_out)
    frustrated_in.loc['Fvalue',:] = np.sum(frustrated_in)
    
    frustrated_out.T.to_csv(saveFilename+'_frustrated_out.csv')
    frustrated_in.T.to_csv(saveFilename+'_frustrated_in.csv')    
    
    
    
    # normalization   
    deltaF_out = pd.DataFrame(frustrated_out.values - frustrated_out.loc[:,str(sorted(ctrl.items()))].values.reshape(-1,1), index = frustrated_out.index, columns = frustrated_out.columns)
    deltaF_in = pd.DataFrame(frustrated_in.values - frustrated_in.loc[:,str(sorted(ctrl.items()))].values.reshape(-1,1), index = frustrated_in.index, columns = frustrated_in.columns)
    
    ind = pd.DataFrame.from_dict(dict(graph.in_degree), orient='index', columns = ['doubleIndegree'])*2
    ind.loc['Fvalue',:] = np.sum(ind).values
    deltaF_in_ = pd.concat([ind, deltaF_in],axis=1,sort=True)
    deltaF_norm_in = pd.concat([ind, pd.DataFrame(deltaF_in_.iloc[:,1:].values/deltaF_in_.iloc[:,0].values.reshape(-1,1), index = deltaF_in_.index, columns = deltaF_in.columns)],axis=1,sort=True)
    deltaF_norm_in.T.to_csv(saveFilename+'_deltaF_normalized_in.csv')

    
    out = pd.DataFrame.from_dict(dict(graph.out_degree), orient='index', columns = ['doubleOutdegree'])*2
    out.loc['Fvalue',:] = np.sum(out)
    deltaF_out_ = pd.concat([out, deltaF_out],axis=1,sort=True)
    deltaF_norm_out = pd.concat([out, pd.DataFrame(deltaF_out_.iloc[:,1:].values/deltaF_out_.iloc[:,0].values.reshape(-1,1), index = deltaF_out_.index, columns = deltaF_out.columns)],axis=1,sort=True)
    deltaF_norm_out.T.to_csv(saveFilename+'_deltaF_normalized_out.csv')


