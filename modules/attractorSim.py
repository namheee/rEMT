# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:25:04 2020

@ author: Namhee Kim
@ We used a Python package ‘PyBoolNet’ for attractor simulation (https://github.com/hklarner/pyboolnet)

"""



import os
import time
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
from typing import List, Dict
import random

from pyboolnet.prime_implicants import create_constants
from pyboolnet.state_transition_graphs import primes2stg

my_env = os.environ.copy()
my_env["PATH"]
Vector1 = Dict[str, str]
Vector2 = List[str]   


#===========================================================
# Random sampling
#===========================================================
def rand_initial_states(num_of_state, num_of_nodes):
    """
    # generate random Boolean initial states
    
    Parameters
    ----------
    num_of_state : int
    num_of_nodes : int     
    """
    if (num_of_nodes <= 32) | (num_of_state*10000 > (2**num_of_nodes)) :     
        rand_int = random.sample(range(2**num_of_nodes), num_of_state)
    else :
        range_num = num_of_state*100000
        rand_int = random.sample(range(range_num), num_of_state)
        print("WARNING : Out of memory!")
        
    s = [bin(x)[2:] for x in list(rand_int)]
    initset = [('0'*(num_of_nodes-len(x))+x) for x in list(s)]   

    return initset

#===========================================================
# Basins of attraction
#===========================================================
def compute_attractor_from_primes(primes, update_mode, initState):
    """
    # compute basin sizes from stg computed by (@pyboolnet)
    # using pyboolnet.state_transition_graphs.primes2stg
    
    Parameters
    ----------
    primes : primes (@pyboolnet)
    update_mode : str ('synchronous' or 'asynchronous')
    initState : List (random sampling)
    """
    stg = primes2stg(primes, update_mode, initState)
    
    attrs_fromSTG = defaultdict()
    for idx, att in enumerate(nx.simple_cycles(stg)):
        attrs_fromSTG[idx] = defaultdict()
        attrs_fromSTG[idx]['attractors'] = tuple(str(x) for x in att)
        attrs_fromSTG[idx]['basinsizes'] = 0
        
        for node in stg.nodes:
            try:
                nx.shortest_path(stg, node, att[0])
                attrs_fromSTG[idx]['basinsizes'] += 1
            except: # there is no path
                continue
            
        attrs_fromSTG[idx]['perc'] = attrs_fromSTG[idx]['basinsizes'] / len(stg.nodes())
    return attrs_fromSTG


#===========================================================
# Phenotype
#===========================================================
def define_phenotype_from_att(primes, phenotype:Vector1, phenotypeAnnot:Vector1, att):
    nodeList = list(primes.keys())
    att_mean = np.mean(att,axis=0)
    annot = 10
    for name, marker in phenotype.items():
        if np.prod([np.rint(np.nextafter((att_mean[nodeList.index(m)]),(att_mean[nodeList.index(m)])+1)) == bool(marker[m]) for m in marker.keys()]) == 1:
            annot = phenotypeAnnot[name]
    return annot

def str2array(x):
    if type(x) == np.str:
        array = np.array([[int(y) for y in x]])
    else:
        array = np.zeros((len(x),len(x[0])))
        for idx, x0 in enumerate(x):
            array[idx,:] = np.array([int(x) for x in x0])
    return array
    
def compute_phenotype(primes, attrs_fromSTG, phenotype, phenotypeAnnot):
    """
    # determine phenotype of an attractor 

    Parameters
    ----------    
    primes : primes (@pyboolnet)
    attrs_fromSTG : the output of 'compute_attractor_from_primes'
    phenotype : Dict ({phenotype:phenotype markers})
    phenotypeAnnot : Dict (simple annotation of phenotype)
    """    
    for idx in range(len(attrs_fromSTG)):
        att = str2array(attrs_fromSTG[idx]['attractors'])
        attrs_fromSTG[idx]['phenotype'] = define_phenotype_from_att(primes, phenotype, phenotypeAnnot, att)

    return attrs_fromSTG



#===========================================================
# Attractor simulation
#===========================================================
def makePhenotypeDF(attrs_dict, phenotypeAnnot):
    df = pd.DataFrame({'phenotype':[[pheno for pheno, phenoIdx in phenotypeAnnot.items() if phenoIdx == x][0] for x in [x['phenotype'] for x in attrs_dict.values()]],
                                     'Ratio':[x['perc'] for x in attrs_dict.values()]})
    pheno_df = df.groupby('phenotype').sum()
    return(pheno_df)

def Simulation(constantDict, primes, update_mode, initState, phenotype, phenotypeAnnot):
    """
    # compute average node activites based on basin sizes
    
    Parameters
    ----------
    constantDict : Dict
    primes : primes (@pyboolnet)
    update_mode : str ('synchronous' or 'asynchronous')
    initState : List (random sampling)
    phenotype : Dict ({phenotype:phenotype markers})
    phenotypeAnnot : Dict (simple annotation of phenotype)
    """
    start = time.time()
    primes_new = create_constants(primes, constantDict,  copy=True)
    attrs_dict = compute_attractor_from_primes(primes_new, update_mode, initState)
    attrs_dict = compute_phenotype(primes_new, attrs_dict, phenotype, phenotypeAnnot)
    print('Attractor simulation time :', time.time()-start)
    pheno_df = makePhenotypeDF(attrs_dict, phenotypeAnnot)
    print(makePhenotypeDF(attrs_dict, phenotypeAnnot))
    
    att_all = []
    for idx, value in enumerate(attrs_dict.values()):
        state01 = value['attractors']
        state_np = np.array([[int(x) for x in state] for state in state01])
        att1_mean = np.mean(state_np, axis=0) * value['perc']
        att_all.append(att1_mean)
    att_ave_pd = pd.DataFrame(np.sum(np.array(att_all),axis=0), index=list(primes.keys()))
    return primes_new, pheno_df, att_ave_pd, attrs_dict

