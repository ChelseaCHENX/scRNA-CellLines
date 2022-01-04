from utils import *
from scipy import linalg
import itertools
from sklearn.mixture import GaussianMixture

def grep_genes(chrom:str, start_pos:int, stop_pos:int,
               gene_pos:pd.DataFrame):
    tmp = gene_pos[gene_pos['chrom'] == chrom]
    tmp = tmp[(tmp['txStart']>start_pos) & (tmp['txEnd']<stop_pos)].drop_duplicates(subset='symbol')
    genes = tmp['symbol'].values

    return genes # np.array

def grep_chrArm(gene:str, 
                gene_pos:pd.DataFrame, chr_pos:pd.DataFrame):   
    chrom, start, end = gene_pos[gene_pos['symbol']==gene][['chrom','txStart','txEnd']].values.flatten()
    idx = [x for x in chr_pos.index if str(chr_pos.loc[x,'chrom'])==chrom and start > chr_pos.loc[x,'chromStart'] and end < chr_pos.loc[x,'chromEnd']]
    if len(idx) == 1:
        idx = idx[0]
    else:
        idx = None
    return idx

def pred_model(X:np.array): # n cells X 1 feature

    # find optimal number of clusters
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, 7)
    # cv_types = ['spherical', 'tied', 'diag', 'full']
    cv_types = ['full']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = GaussianMixture(n_components=n_components,
                                        covariance_type=cv_type)
            gmm.fit(X)
            bic.append(gmm.bic(X))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
    return best_gmm

def calcMovingAvg(cts:pd.DataFrame, gene_pos:pd.DataFrame, win=100, inputForm='rawcounts'):
    '''
    cts: gene X cells: raw counts, not Log or Normed
    win: number of genes to average in one window
    '''
    
    order = ['chr%d'%(i) for i in range(1,23)] + ['chrX', 'chrY']
    for n, chrom in enumerate(order):
        tmp = gene_pos[gene_pos['chrom']==chrom].sort_values(by=['txStart'])
        if n == 0:
            res = tmp
        elif n > 0:
            res = pd.concat((res, tmp), axis=0)

    cts = cts.loc[[g for g in res['symbol'] if g in cts.index],:]
    if inputForm == 'rawcounts':
        X = cts.values 
        E = np.log2(1+1e4*X/np.sum(X, axis=0))
    elif inputForm == 'log':
        E = cts.values 
    ER = E - np.mean(E, axis=1, keepdims=True) # each gene among cells, sum=0   
    ER_MvAvg = np.apply_along_axis(lambda x:np.convolve(x, np.ones((win,))/win, mode='same'), 0, ER)

    genes_s, cells = list(cts.index), cts.columns
    return genes_s, cells, ER_MvAvg

    

if __name__ == "__main__":

    wkdir = '/bgfs/alee/chelsea/projects/10X/AL1/codes'
    os.chdir(wkdir)

    # Prep files
    chr_pos = pd.read_csv('../data/hg38_chromPos.csv', index_col=0)
    gene_pos = pd.read_csv('../data/CNV/gene_pos.csv')

    adata = sc.read('../data/scVelo/AL1_ref_raw.h5ad')
    CellTypes = ['MCF7','T47D WT','T47D KO',  'MM134','SUM44','BCK4', 'MCF10A','HEK293'] 
    CNVDict = {}

    min_genes_per_arm = 100
    min_cl_cells = 20
    min_cl_prob = 0.95

    for cell in CellTypes:
        CNVDict[cell] = {}
        cnv = None
        n = -1
        chrs = []
        tmp = adata_ref_raw[adata_ref_raw.obs['CellType']==cell]
        
        X = np.array(tmp.X.todense()).T # genes X cell
        cts_raw = pd.DataFrame(X, index=tmp.var.index, columns=tmp.obs.index) # gene X cell
        
        genes, cells, ER_MvAvg = calcMovingAvg(cts_raw, gene_pos, win=100) # gene X cell
        cts = pd.DataFrame(ER_MvAvg, index=genes, columns=cells)
        CNVDict[cell]['WindowCts'] = cts.T # cell X genes
        
        CNVDict[cell]['gmm'] = {}
        for arm in chr_pos.index:
            
            chrom, start_pos, stop_pos = chr_pos.loc[arm,:] 
            genes = grep_genes(chrom,start_pos,stop_pos,gene_pos)
            genes_s = [g for g in genes if g in cts.index]
            
            if len(genes_s) < min_genes_per_arm:
                print('Low cov: {} {}'.format(cell, chrom))
            if len(genes_s) > min_genes_per_arm:
                x = np.mean(cts.loc[genes_s,:], axis=0).values 
                X = np.reshape(x, (len(x),1))            
                best_gmm = pred_model(X)
                
                if best_gmm.n_components > 1:
                    predicts = best_gmm.predict(X)
                    probs = best_gmm.predict_proba(X)
                    
                    for i in range(best_gmm.n_components):
                        cells_cl = np.where(predicts==i)[0]
                        probs_cl = probs[cells_cl,i].flatten()
                        if sum(probs_cl > min_cl_prob) < min_cl_cells:  
                            break 
                        elif sum(probs_cl > min_cl_prob) > min_cl_cells: 
                            if i < best_gmm.n_components-1:
                                continue
                            if i==(best_gmm.n_components-1):  
                                print('\n',cell, arm, len(genes_s), best_gmm.n_components, Counter(best_gmm.predict(X)))
                                n += 1 
                                chrs.append(arm)
                                if n == 0:
                                    cnv = x
                                elif n > 0:
                                    cnv = np.vstack((cnv, x))
                                CNVDict[cell]['gmm'][arm] = best_gmm

        if cnv is None or len(cnv.shape) < 2:
            continue

        cnv_df = pd.DataFrame(cnv, index=chrs, columns=cts.columns)  # Arm X cells
        CNVDict[cell]['ArmCts'] = cnv_df.T # cell X Arms
        
        # heatmap (cell X Arm)
        width = len(chrs)*0.7
        g = sns.clustermap(cnv_df.T, method='ward', metric='euclidean', 
                        cmap='coolwarm', figsize=(width, 5))
        
        g.ax_heatmap.set_yticklabels('')
        
        g.ax_heatmap.set_xlabel(cell)
        xticklabels = cnv_df.T.columns[np.array(g.dendrogram_col.reordered_ind)]
        g.ax_heatmap.set_xticks([x+0.5 for x in range(len(xticklabels))])
        g.ax_heatmap.set_xticklabels(xticklabels, rotation=60)

        cbar = g.cax
        x0,y0,x1,y1 = np.array(g.ax_heatmap.get_position()).flatten()
        cbar.set_position([1.3, y0, 0.2, y1-y0])

        CNVDict[cell]['ArmCtsHeatmap'] = g 
    
    