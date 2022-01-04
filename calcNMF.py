import numpy as np
from sklearn.decomposition import NMF 

def calcNMF(E:np.array, n_genes=50, kmin=5, kmax=15):
    '''
    E: cell X genes: log2(cts+1)
    Relative Expression Mat = tmp => each gene: centerized
    calc NMF score `http://garmiregroup.org/nmfem`
    return genes with `num` top NMF score
    '''       
    Er = E - np.mean(E, axis=0)
    
    Er = np.clip(Er, a_min=0, a_max=None)
    Const = np.full((Er.shape[0],1), np.mean(Er, axis=None))
    Er = np.hstack((Er, Const))

    gOrder_list = []

    for k in range(kmin, kmax):
        model = NMF(n_components=k, init='random', random_state=0, alpha=0)
        W = model.fit_transform(Er)
        H = model.components_ # k comps X N genes
        H = np.divide(H, H[:,[-1]]) # normalize using components of the const. col
        D_score = np.apply_along_axis(lambda x: max(x)-min(x), 0, H)[:-1]
        gOrder = np.argsort(D_score)[::-1]
        gOrder_list.append(gOrder[0:n_genes])
        
    return np.array(gOrder_list) 