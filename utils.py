import scvelo as scv
import scanpy as sc
# import scrublet as scr

import pickle
import anndata 
# from MulticoreTSNE import MulticoreTSNE as TSNE

import pandas as pd
from pandas.api.types import CategoricalDtype
import numpy as np
from math import *
from collections import Counter

from scipy.stats import zscore
from scipy.spatial.distance import cdist
from scipy.linalg import svd
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA

import subprocess
from os.path import isfile
import os
import sys
import re
import glob
import logging as logg

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

def run(cmd):
    proc = subprocess.Popen([cmd], shell=True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
    )
    stdout, stderr = proc.communicate()
    return proc.returncode, stdout, stderr

def getarg(argv, index, default=None):
    if len(argv) > index:
        ret = argv[index]
    else:
        ret = default
    return ret

def addGene2Obs(genes, adata_ref, adata_query, layer=None, prefix='Counts'):
    '''
    gene: `list` of strings
    adata_ref: raw adata with all genes detected
    adata_query: adata to be added value to plot 
    layer: layer of adata_ref to add. 'ambiguous', 'matrix', 'spliced', 'unspliced', 'Ms', 'Mu'. Default: adata.X, i.e., spliced counts'''
    if isinstance(genes, list):
        cells = adata_query.obs.index.tolist()
        genes = list(set(genes) & set(adata_ref.var.index))
        genes_s = [gene for gene in genes if '%s_%s'%(prefix, gene) not in adata_query.obs.columns]
        
        idx_cells = [np.where(adata_ref.obs.index == cell)[0][0] for cell in cells]
        idx_genes = [np.where(adata_ref.var.index == gene)[0][0] for gene in genes_s]
        print('genes number %d'%(len(genes)))
        if len(idx_genes) > 0:
            expr = pd.DataFrame(adata_ref.X.todense()[idx_cells, :][:, idx_genes], columns=['%s_%s'%(prefix, gene) for gene in genes_s], index=cells)
            nobs = adata_query.obs.join(expr)
            return nobs
        else:
            print('Non obs added')
            return adata_query.obs
    else:
        print('Input genes not a str list')
        return adata_query.obs

    
def gather_degs(adata):
    res = adata.uns['rank_genes_groups']
    cats = adata.obs[res['params']['groupby']].dtype.categories
    print(cats)
    deg_dict = {}

    for cat in cats:
        for n, name in enumerate(['scores','logfoldchanges','pvals','pvals_adj']):
            tmp = res[name][cat]
            if n == 0:
                out = tmp
            elif n > 0:
                out = np.vstack((out, tmp))

        outdf = pd.DataFrame(out, index=['scores','logfoldchanges','pvals','pvals_adj'], columns=res['names'][cat]).T
        deg_dict[cat] = outdf        
    return deg_dict

def corr2_coeff(A,B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Finally get corr coeff
    return np.dot(A_mA,B_mB.T)/np.sqrt(np.dot(ssA[:,None],ssB[None]))    

def gather_degs(adata):
    res = adata.uns['rank_genes_groups']
    cats = adata.obs[res['params']['groupby']].dtype.categories
    print(cats)
    deg_dict = {}

    for cat in cats:
        for n, name in enumerate(['scores','logfoldchanges','pvals','pvals_adj']):
            tmp = res[name][cat]
            if n == 0:
                out = tmp
            elif n > 0:
                out = np.vstack((out, tmp))

        outdf = pd.DataFrame(out, index=['scores','logfoldchanges','pvals','pvals_adj'], columns=res['names'][cat]).T
        deg_dict[cat] = outdf       
    return deg_dict  
    
def read_gmt(fpath):
    gdict = {}
    
    with open(fpath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            items = line.strip().rstrip().split('\t')
            gset = items[0].lower()
            genes = items[2:]
            gdict[gset] = genes            
    return gdict

def write_gmt(gdict, fpath):  
    with open(fpath, 'w') as f:
        for key, val in gdict.items():
            val_string = '\t'.join(val)
            f.write('%s\thttp\t%s\n'%(key, val_string))            
    print('finished writing %d gene sets to %s'%(len(gdict), fpath))


# dissect clustering dendrograms
def get_leaves(k, lk_df):
    '''
    k: cluster number `father node`
    lk_df: linkage df
    return: a list of leaves by original order
    '''
    res = []
    N = lk_df.shape[0]+1
    
    if k > N-1:

        n1, n2 = lk_df[int(k)-N,[0,1]]
        lv1 = get_leaves(n1,lk_df)
        lv2 = get_leaves(n2,lk_df)
        res.extend(lv1)
        res.extend(lv2)
        return res 
    
    elif k < N:
        return [int(k)]

def get_clade_labels(labels, lk_df, depth=1):
    '''
    labels: list, labels of the original nodes
    depth `h`:int, h>0, the last h step in linkage merging, corresponding to (h+1) clusters
    lk_df: linkage df
    '''
    N = lk_df.shape[0]+1
    res_dict = {}
    
    idx = 2*N-2-depth
    nodes = lk_df[(N-1-depth): N,[0,1]].flatten()
    nodes = [n for n in nodes if n < 2*N-depth-1]
        
#     res_dict = {n: [labels[int(x)] for x in [get_leaves(n, lk_df)]] for n in nodes}
    res_dict = {int(n): [labels[int(x)] for x in get_leaves(n, lk_df)] for n in nodes}
    return res_dict


def calcPairwisePcDist(X:np.array, n_pcs=30, metric='euclidean'):
    '''
    X: obs x features, `m x n` np.array; not centered or scaled for each feature
    n_pcs: 
    metric: pairwise distance metric
    return a `m x m` pairwise distance matrix using n_pcs=n_pcs
    '''
    from sklearn.decomposition import PCA
    from scipy.spatial.distance import pdist, squareform

    X = X - np.mean(X, axis=0) # obs x features
    pca = PCA(n_components=n_pcs)
    pca.fit(X)

    Y = np.dot(pca.components_, X.T).T # obs x features
    D = squareform(pdist(Y, metric=metric))

    return D 