# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:46:54 2022

@author: Namhee Kim
"""

from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd

def phenotypeAnnot_diff(x):
    if x>=0.5: return (-1)
    elif x<(-0.5): return (1)
    else: return (0)

def perturb_analysis(perturb_p, perturb_s, save_perturbname):
    perturb_p.loc[:,'Desired1'] = [1,1,1,0,0,0,1,0,0] #rEMT
    perturb_p.loc[:,'Desired2'] = [1,1,1,0,0,0,0,0,0] #chemosensitive rEMT
        
    dT1 = perturb_p.loc[:,sorted(set(perturb_p.columns)-set(['Desired2']))]
    didx1 = list(dT1.columns).index('Desired1')
    dc_T1_ = cosine_similarity(dT1.T)
    dT1.loc['cosineSim1',:] = dc_T1_[didx1,:]
    
    dT2 = perturb_p.loc[:,sorted(set(perturb_p.columns)-set(['Desired1']))]
    didx2 = list(dT2.columns).index('Desired2')
    dc_T2_ = cosine_similarity(dT2.T)
    dT2.loc['cosineSim2',:] = dc_T2_[didx2,:]
    
    dT3 = pd.concat([dT1.loc['cosineSim1',:], dT2.loc['cosineSim2',:]], axis=1, sort=True).T
    perturb_p.loc['cosineSim1',:] = dT3.loc['cosineSim1',perturb_p.columns]
    perturb_p.loc['cosineSim2',:] = dT3.loc['cosineSim2',perturb_p.columns]
    perturb_p = perturb_p.fillna(0)
    perturb_p.loc['Diff',:] =  perturb_p.loc['Ecadherin',:] -  perturb_p.loc['ZEB1',:]
    perturb_p.loc['Diff_pheno',:] = [phenotypeAnnot_diff(x) for x in (perturb_p.loc['Ecadherin',:] -  perturb_p.loc['ZEB1',:]).values]    
    
    ## rEMT score
    perturb_p.loc['weight',:] = (perturb_p.loc['Diff',:] + 1)/2
    perturb_p.loc['rEMTscore',:] = perturb_p.loc['cosineSim1',:].values * perturb_p.loc['weight',:].values
    perturb_p.loc['chemo_rEMTscore',:] = perturb_p.loc['cosineSim2',:].values * perturb_p.loc['weight',:].values
    
    # save result
    perturb_result = pd.concat([perturb_p, perturb_s], sort=True).T        
    perturb_result.to_csv(save_perturbname)    