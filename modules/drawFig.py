# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 10:44:23 2022

@author: Namhee Kim

# To determine a best-fit surface, we used numpy and scipy packages according to the reference.
# Please refer to (Ref. https://gist.github.com/amroamroamro).

"""

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

def removeN1(x): # remove phenotype markers and mutation
    for s in ['ZEB1','Ecadherin','RAS']:
        if s in x:
            return False
            break
        else: continue
    return True

def draw_surface3D(df, axis_name, surface_color):
    data = np.array(df.loc[:,axis_name])
    
    min_ = np.min(data, axis=0)
    max_ = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(min_[0], max_[0], 20), np.linspace(min_[1], max_[1], 20))
    XX = X.flatten()
    YY = Y.flatten()
        
    # best-fit quadratic curve (2nd-order)
    A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
        
    # evaluate it on a grid
    Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
        
    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z-0.01, rstride=2, cstride=2, alpha=0.2, cmap=surface_color)

    palette = {'#9E0142':'#9E0142','#F39001':'#F39001','#5E4FA2':'#5E4FA2',"#98FB98": "#98FB98", "#3CB371":'#3CB371'}
    for color in set(df.Pheno):    
        dataset = df.loc[(df.Pheno==color).values,axis_name].values
        print(dataset.shape)
        ax.scatter(dataset[:,0], dataset[:,1], dataset[:,2], c=palette[color], alpha=0.5, s=150)
        ax.view_init(15, 60)

    plt.xlabel(axis_name[0])
    plt.ylabel(axis_name[1])
    ax.set_zlabel(axis_name[2])
