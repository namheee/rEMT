# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:49:03 2022

@author: Namhee Kim
"""

import numpy as np
from collections import defaultdict

def compute_networkStability(attrs_dict, graph, nodeList):
    """
    Parameters
    ----------
    attrs_dict : Dict
    graph : nx.Digraph
    nodeList : List

    Returns
    -------
    frustration_dict : Dict
        frustration of each network state.
    network_stability : Dict
        summarize the network stability; attractor stability and frustration of a major attractor.

    """
    network_stability = {'att_stability':0,'F_major_att':0,'major_att_size':0}
    major_size = np.max([value['perc'] for value in attrs_dict.values()])
    major_atti = [key for key, value in attrs_dict.items() if value['perc']==major_size][0]
    frustration_dict = defaultdict(dict)
    
    for idx, value in enumerate(attrs_dict.values()):
        # 1. attractor stability
        # ref. Rachdi, Mustapha, et al. "Entropy as a robustness marker in genetic regulatory networks." Entropy 22.3 (2020): 260.
        network_stability['att_stability'] += -1*attrs_dict[idx]['perc']*np.log(attrs_dict[idx]['perc'])
        
        # 2. frustration of network state
        # ref. Tripathi, Shubham, David A. Kessler, and Herbert Levine. "Biological networks regulating cell fate choice are minimally frustrated." Physical Review Letters 125.8 (2020): 088101.
        state01 = value['attractors']
        frustration_dict[idx]['frustration'] = 0
        for state in state01:
            state_np = [-1 if int(x) == 0 else int(x) for x in state]
            frustration_dict[idx]['perc'] = attrs_dict[idx]['perc']
            frustration_dict[idx]['frustration'] += np.sum([-1*list(graph.adj[x][y]['sign'])[0]* state_np[nodeList.index(x)]* state_np[nodeList.index(y)] for x,y in graph.edges()])

        frustration_dict[idx]['frustration'] = frustration_dict[idx]['frustration'] / len(state01) #cyclic attractor
    network_stability['F_major_att'] = frustration_dict[major_atti]['frustration']
    network_stability['major_att_size'] = frustration_dict[major_atti]['perc']
    return frustration_dict, network_stability


    